module Ginzburg_Landau_Base
   ! Standard Library.
   use stdlib_stats_distribution_normal, only: normal => rvs_normal
   use stdlib_optval, only : optval
   use stdlib_math, only : linspace
   use stdlib_linalg, only : eye
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
   character(len=128), parameter :: this_module = 'Ginzburg_Landau_Base'
   
   public  :: L, nx, N, dx
   public  :: nu, gamma, mu_0, c_mu, mu_2, mu
   public  :: rk_b, x_b, s_b, rk_c, x_c, s_c
   public  :: B, CT, weight, weight_mat
   public  :: BBTW, CTCW
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
   real(wp)               :: mu(nx)

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
   real(wp)               :: weight_mat(N,N)    ! integration weights matrix
   real(wp)               :: weight_flat(N**2)    ! integration weights flat
   real(wp)               :: BBTW(N,N)
   real(wp)               :: CTCW(N,N)
   ! Data matrices for Riccati
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
      real(wp), dimension(2*nx) :: mean, std
      normalize = optval(ifnorm,.true.)
      mean = 0.0_wp
      std  = 1.0_wp
      self%state = normal(mean,std)
      if (normalize) then
         alpha = self%norm()
         call self%scal(1.0/alpha)
      endif
      return
   end subroutine rand

   !------------------------------------------------------
   !-----     TYPE BOUND PROCEDURES FOR LR STATES    -----
   !------------------------------------------------------

   subroutine initialize_LR_state(self, U, S, rk, rkmax, if_rank_adaptive)
      class(LR_state),            intent(inout) :: self
      class(abstract_vector_rdp), intent(in)    :: U(:)
      real(wp),                   intent(in)    :: S(:,:)
      integer,                    intent(in)    :: rk
      integer, optional,          intent(in)    :: rkmax
      logical, optional,          intent(in)    :: if_rank_adaptive
      logical                                   :: ifrk

      ! internals
      class(abstract_vector_rdp), allocatable   :: Utmp(:)
      real(wp), allocatable :: R(:, :)
      integer :: i, m, rka, info
      character(len=128) :: msg

      ifrk = optval(if_rank_adaptive, .false.)

      select type (U)
      type is (state_vector)
         m = size(U)
         call assert_shape(S, [m,m], 'S', this_module, 'initialize_LR_state')
         ! optional size argument
         if (present(rkmax)) then
            if (rkmax < rk) then
               call stop_error('rkmax < rk!', this_module, 'initialize_LR_state')
            end if
            self%rk = rk
            rka = rkmax
            if (ifrk) then
               if (rkmax==rk) then
                  call stop_error('rkmax must be larger than rk for rank-adaptive DLRA!', this_module, 'initialize_LR_state')
               end if
               write(msg,'(A,I0,A)') 'Rank-adaptivity enabled. Computation will begin with X%rk = ', self%rk+1, '.'
               call logger%log_information(msg, module=this_module, procedure='initialize_LR_state')
            end if
         else
            self%rk = rk
            if (ifrk) then
               rka = rk + 1
            else
               rka = rk
            end if
         end if

         ! allocate & initialize
         allocate(self%U(rka), source=U(1)); call zero_basis(self%U)
         allocate(self%S(rka,rka)); self%S = 0.0_wp
         write(msg,'(3(A,I0),A)') 'size(X%U) = [ ', rka,' ], X%rk = ', self%rk, ', size(U0) = [ ', m,' ]'
         call logger%log_information(msg, module=this_module, procedure='initialize_LR_state')
         ! copy inputs
         if (self%rk > m) then   ! copy the full IC into self%U
            call copy_basis(self%U(:m), U)
            self%S(:m,:m) = S
            write(msg,'(4X,A,I0,A)') 'Transfer the first ', m, ' columns of U0 to X%U.'
            call logger%log_information(msg, module=this_module, procedure='initialize_LR_state')
         else  ! fill the first self%rk columns of self%U with the first self%rk columns of the IC
            call copy_basis(self%U(:self%rk), U(:self%rk))
            self%S(:self%rk,:self%rk) = S(:self%rk,:self%rk)
            write(msg,'(4X,A,I0,A)') 'Transfer all ', m, ' columns of U0 to X%U.'
            call logger%log_information(msg, module=this_module, procedure='initialize_LR_state')
         end if
         ! top up basis (to rka for rank-adaptivity) with orthonormal columns if needed
         if (rka > m) then
            write(msg,'(4X,A,I0,A)') 'Fill remaining ', rka-m, ' columns with orthonormal noise orthonormal to U0.'
            call logger%log_information(msg, module=this_module, procedure='initialize_LR_state')
            allocate(Utmp(rka-m), source=U(1))
            do i = 1, rka-m
               call Utmp(i)%rand()
            end do
            allocate(R(rka-m,rka-m)); R = 0.0_wp
            call orthogonalize_against_basis(Utmp, self%U, info)
            call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='initialize_LR_state')
            call qr(Utmp, R, info)
            call check_info(info, 'qr', module=this_module, procedure='initialize_LR_state')
            call copy_basis(self%U(m+1:), Utmp)
         end if
      end select
      return
   end subroutine initialize_LR_state

end module Ginzburg_Landau_Base
