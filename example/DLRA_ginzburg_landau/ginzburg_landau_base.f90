module Ginzburg_Landau_Base
   ! Standard Library.
   use stdlib_optval, only : optval
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov, only: wp => dp
   use LightKrylov_Logger
   use LightKrylov_Utils, only : assert_shape
   use LightKrylov_AbstractVectors
   ! LightROM
   use LightROM_AbstractLTIsystems ! LR_state
   implicit none

   private :: this_module
   character*128, parameter :: this_module = 'Ginzburg_Landau_Base'
   
   public  :: L, nx, dx
   public  :: nu, gamma, mu_0, c_mu, mu_2, mu
   public  :: rk_b, x_b, s_b, rk_c, x_c, s_c
   public  :: B, CT, weight, weight_mat
   public  :: N, BBTW_flat, CTCW_flat
   public  :: Qc, Rinv, CTQcCW_mat, BRinvBTW_mat

   
   !-------------------------------
   !-----     PARAMETERS 1    -----
   !-------------------------------

   ! Mesh related parameters.
   real(wp), parameter :: L  = 50.0_wp ! Domain length
   integer,  parameter :: nx = 128     ! Number of grid points (excluding boundaries).
   real(wp)            :: dx           ! Grid size.

   !-------------------------------------------
   !-----     LIGHTKRYLOV VECTOR TYPE     -----
   !-------------------------------------------

   type, extends(abstract_vector_rdp), public :: state_vector
      real(wp) :: state(2*nx) = 0.0_wp
   contains
      private
      procedure, pass(self), public :: zero
      procedure, pass(self), public :: dot
      procedure, pass(self), public :: scal
      procedure, pass(self), public :: axpby
      procedure, pass(self), public :: rand
      procedure, pass(self), public :: get_size
   end type state_vector

   !-------------------------------------------------------
   !-----     LIGHTKRYLOV SYM LOW RANK STATE TYPE     -----
   !-------------------------------------------------------

   type, extends(abstract_sym_low_rank_state_rdp), public :: LR_state
   contains
      private
      procedure, pass(self), public :: initialize_LR_state
   end type LR_state

   !-------------------------------
   !-----     PARAMETERS 2    -----
   !-------------------------------

   ! Physical parameters.
   complex(wp), parameter :: nu    = cmplx(2.0_wp, 0.2_wp, wp)
   complex(wp), parameter :: gamma = cmplx(1.0_wp, -1.0_wp, wp)
   real(wp),    parameter :: mu_0  = 0.38_wp
   real(wp),    parameter :: c_mu  = 0.2_wp
   real(wp),    parameter :: mu_2  = -0.01_wp
   real(wp)               :: mu(1:nx)

   ! Input-Output system parameters
   real(wp)               :: weight(2*nx)       ! integration weights
   integer,  parameter    :: rk_b = 2           ! number of inputs to the system
   real(wp), parameter    :: x_b = -11.0_wp     ! location of input Gaussian
   real(wp), parameter    :: s_b = 1.0_wp       ! variance of input Gaussian
   type(state_vector)     :: B(rk_b)
   real(wp), parameter    :: x_c = sqrt(-2.0_wp*(mu_0 - c_mu**2)/mu_2) ! location of input Gaussian
   real(wp), parameter    :: s_c = 1.0_wp       ! variance of input Gaussian
   integer,  parameter    :: rk_c = 2           ! number of outputs to the system
   type(state_vector)     :: CT(rk_c)
   real(wp)               :: Qc(rk_c,rk_c)
   real(wp)               :: Rinv(rk_b,rk_b)

   ! Data matrices for RK lyap
   integer,  parameter    :: N = 2*nx           ! Number of grid points (excluding boundaries).
   real(wp)               :: weight_mat(N**2)   ! integration weights
   real(wp)               :: BBTW_flat(N**2)
   real(wp)               :: CTCW_flat(N**2)
   ! Data matrices for Riccatis
   real(wp)               :: CTQcCW_mat(N,N)
   real(wp)               :: BRinvBTW_mat(N,N)

contains

   !=========================================================
   !=========================================================
   !=====                                               =====
   !=====     LIGHTKRYLOV MANDATORY IMPLEMENTATIONS     =====
   !=====                                               =====
   !=========================================================
   !=========================================================

   !----------------------------------------------------
   !-----     TYPE-BOUND PROCEDURE FOR VECTORS     -----
   !----------------------------------------------------

   subroutine zero(self)
      class(state_vector), intent(inout) :: self
      self%state = 0.0_wp
      return
   end subroutine zero

   real(wp) function dot(self, vec) result(alpha)
      ! weighted inner product
      class(state_vector),        intent(in) :: self
      class(abstract_vector_rdp), intent(in) :: vec
      select type(vec)
      type is(state_vector)
         alpha = dot_product(self%state, weight*vec%state)
      end select
      return
   end function dot

   integer function get_size(self) result(N)
     class(state_vector), intent(in) :: self
     N = 2*nx
     return
   end function get_size

   subroutine scal(self, alpha)
      class(state_vector), intent(inout) :: self
      real(wp),            intent(in)    :: alpha
      self%state = self%state * alpha
      return
   end subroutine scal

   subroutine axpby(self, alpha, vec, beta)
      class(state_vector),        intent(inout) :: self
      class(abstract_vector_rdp), intent(in)    :: vec
      real(wp),                   intent(in)    :: alpha, beta
      select type(vec)
      type is(state_vector)
         self%state = alpha*self%state + beta*vec%state
      end select
      return
   end subroutine axpby

   subroutine rand(self, ifnorm)
      class(state_vector), intent(inout) :: self
      logical, optional,   intent(in)    :: ifnorm
      ! internals
      logical :: normalize
      real(wp) :: alpha
      normalize = optval(ifnorm,.true.)
      call random_number(self%state)
      if (normalize) then
         alpha = self%norm()
         call self%scal(1.0/alpha)
      endif
      return
   end subroutine rand

   !------------------------------------------------------
   !-----     TYPE BOUND PROCEDURES FOR LR STATES    -----
   !------------------------------------------------------

   subroutine initialize_LR_state(self, U, S, rk, rkmax)
      class(LR_state),            intent(inout) :: self
      class(abstract_vector_rdp), intent(in)    :: U(:)
      real(wp),                   intent(in)    :: S(:,:)
      integer,                    intent(in)    :: rk
      integer, optional,          intent(in)    :: rkmax

      ! internals
      real(wp), allocatable :: R(:, :)
      integer :: i, n, rka, info

      n = size(U)
      call assert_shape(S, [n,n], "initialize_LR_state", "S")

      ! optional size argument
      if (present(rkmax)) then
         self%rk = rk - 1
         rka = rkmax
      else
         self%rk = rk
         rka = rk + 1
      end if

      select type (U)
      type is (state_vector)
         ! allocate & initialize
         allocate(self%U(rka), source=U(1)); call zero_basis(self%U)
         allocate(self%S(rka,rka)); self%S = 0.0_wp
         ! copy inputs
         if (self%rk > n) then   ! copy the full IC into self%U
            call copy_basis(self%U(1:n), U)
            self%S(1:n,1:n) = S
         else  ! fill the first self%rk columns of self%U with the first self%rk columns of the IC
            call copy_basis(self%U(1:self%rk), U(1:self%rk))
            self%S(1:self%rk,1:self%rk) = S(1:self%rk,1:self%rk)
         end if
         ! top up basis (to rka for rank-adaptivity) with orthonormal columns if needed
         if (rka > n) then
            do i = n+1, rka
               call self%U(i)%rand()
            end do
            allocate(R(rka,rka)); R = 0.0_wp
            call qr(self%U, R, info)
            call check_info(info, 'qr', module=this_module, procedure='initialize_LR_state')
         end if
      end select
      return
   end subroutine initialize_LR_state

end module Ginzburg_Landau_Base
