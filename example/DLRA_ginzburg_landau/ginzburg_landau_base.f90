module Ginzburg_Landau_Base
   ! Standard Library.
   use stdlib_stats_distribution_normal, only: normal => rvs_normal
   use stdlib_optval, only : optval
   use stdlib_math, only : linspace
   use stdlib_io_npy, only: save_npy
   use stdlib_strings, only: replace_all
   use stdlib_linalg, only: svdvals, eye
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov_Logger
   use LightKrylov_Utils, only : assert_shape
   use LightKrylov_Constants, only : one_rdp, zero_rdp
   use LightKrylov_AbstractVectors
   ! LightROM
   use LightROM_AbstractLTIsystems ! LR_state
   implicit none

   private :: this_module
   character(len=*), parameter :: this_module = 'Ginzburg_Landau_Base'
   
   public  :: L, nx, N, dx
   public  :: nu, gamma, mu_0, c_mu, mu_2, mu
   public  :: rk_b, x_b, s_b, rk_c, x_c, s_c
   public  :: B, CT
   public  :: weight, weight_mat
   public  :: BBTW, CTCW
   public  :: Qc, Rinv, Qe, Vinv
   public  :: CTQcCW, BRinvBTW, CTVinvCW, BQeBTW
   
   !-------------------------------
   !-----     PARAMETERS 1    -----
   !-------------------------------

   ! Mesh related parameters.
   real(dp), parameter :: L  = 50.0_dp ! Domain length
   integer,  parameter :: nx = 128     ! Number of grid points (excluding boundaries).
   real(dp)            :: dx           ! Grid size.

   !--------------------------------------------
   !-----     LIGHTKRYLOV VECTOR TYPES     -----
   !--------------------------------------------

   type, extends(abstract_vector_rdp), public :: state_vector
      real(dp) :: state(2*nx) = 0.0_dp
   contains
      private
      procedure, pass(self), public :: zero
      procedure, pass(self), public :: dot
      procedure, pass(self), public :: scal
      procedure, pass(self), public :: axpby
      procedure, pass(self), public :: rand
      procedure, pass(self), public :: get_size
   end type state_vector

   type, extends(abstract_vector_rdp), public :: state_error_vector
      type(state_vector) :: state
      type(state_vector) :: error
   contains
      private
      procedure, pass(self), public :: zero => se_zero
      procedure, pass(self), public :: dot => se_dot
      procedure, pass(self), public :: scal => se_scal
      procedure, pass(self), public :: axpby => se_axpby
      procedure, pass(self), public :: rand => se_rand
      procedure, pass(self), public :: get_size => se_get_size
   end type state_error_vector

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
   integer,  parameter    :: N = 2*nx           ! Number of grid points (excluding boundaries).
   complex(dp), parameter :: nu    = cmplx(2.0_dp, 0.2_dp, dp)
   complex(dp), parameter :: gamma = cmplx(1.0_dp, -1.0_dp, dp)
   real(dp),    parameter :: mu_0  = 0.38_dp
   real(dp),    parameter :: c_mu  = 0.2_dp
   real(dp),    parameter :: mu_2  = -0.01_dp
   real(dp)               :: mu(nx)

   ! Input-Output system parameters
   real(dp)               :: weight(N)          ! integration weights
   integer,  parameter    :: rk_b = 2           ! number of inputs to the system
   real(dp), parameter    :: x_b = -11.0_dp     ! location of input Gaussian
   real(dp), parameter    :: s_b = 1.0_dp       ! variance of input Gaussian
   type(state_vector)     :: B(rk_b)
   real(dp), parameter    :: x_c = sqrt(-2.0_dp*(mu_0 - c_mu**2)/mu_2) ! location of input Gaussian
   real(dp), parameter    :: s_c = 1.0_dp       ! variance of input Gaussian
   integer,  parameter    :: rk_c = 2           ! number of outputs to the system
   type(state_vector)     :: CT(rk_c)
   real(dp)               :: Qc(rk_c,rk_c)      ! State energy
   real(dp)               :: Rinv(rk_b,rk_b)    ! Inverse control cost
   real(dp)               :: Vinv(rk_c,rk_c)    ! Variance of the sensor noise
   real(dp)               :: Qe(rk_b,rk_b)      ! Variance of the actuator noise

   ! Data matrices for RK lyap
   real(dp)               :: weight_mat(N,N)    ! integration weights matrix
   real(dp)               :: weight_flat(N**2)  ! integration weights flat
   real(dp)               :: BBTW(N,N)
   real(dp)               :: CTCW(N,N)
   ! Data matrices for Riccati
   real(dp)               :: CTVinvCW(N,N)
   real(dp)               :: BQeBTW(N,N)
   real(dp)               :: BRinvBTW(N,N)
   real(dp)               :: CTQcCW(N,N)

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
      self%state = 0.0_dp
   end subroutine zero

   real(dp) function dot(self, vec) result(alpha)
      ! weighted inner product
      class(state_vector),        intent(in) :: self
      class(abstract_vector_rdp), intent(in) :: vec
      select type(vec)
      type is(state_vector)
         alpha = dot_product(self%state, weight*vec%state)
      class default
         call type_error('vec', 'state_vector', 'IN', this_module, 'dot')
      end select
   end function dot

   integer function get_size(self) result(N)
     class(state_vector), intent(in) :: self
     N = 2*nx
   end function get_size

   subroutine scal(self, alpha)
      class(state_vector), intent(inout) :: self
      real(dp),            intent(in)    :: alpha
      self%state = self%state * alpha
   end subroutine scal

   subroutine axpby(alpha, vec, beta, self)
      class(state_vector),        intent(inout) :: self
      class(abstract_vector_rdp), intent(in)    :: vec
      real(dp),                   intent(in)    :: alpha, beta
      select type(vec)
      type is(state_vector)
         self%state = beta*self%state + alpha*vec%state
      class default
         call type_error('vec', 'state_vector', 'IN', this_module, 'axpby')
      end select
   end subroutine axpby

   subroutine rand(self, ifnorm)
      class(state_vector), intent(inout) :: self
      logical, optional,   intent(in)    :: ifnorm
      ! internals
      logical :: normalize
      real(dp) :: alpha
      real(dp), dimension(2*nx) :: mean, std
      normalize = optval(ifnorm,.true.)
      mean = 0.0_dp
      std  = 1.0_dp
      self%state = normal(mean,std)
      if (normalize) then
         alpha = self%norm()
         call self%scal(1.0_dp/alpha)
      endif
   end subroutine rand

   !----------------------------------------------------
   !-----     TYPE-BOUND PROCEDURE FOR VECTORS     -----
   !----------------------------------------------------

   subroutine se_zero(self)
      class(state_error_vector), intent(inout) :: self
      call self%state%zero()
      call self%error%zero()
   end subroutine se_zero

   real(dp) function se_dot(self, vec) result(alpha)
      ! weighted inner product
      class(state_error_vector),  intent(in) :: self
      class(abstract_vector_rdp), intent(in) :: vec
      select type(vec)
      type is(state_error_vector)
         alpha = self%state%dot(vec%state) + self%error%dot(vec%error)
      class default
         call type_error('vec', 'state_error_vector', 'IN', this_module, 'dot')
      end select
   end function se_dot

   integer function se_get_size(self) result(N)
     class(state_error_vector), intent(in) :: self
     N = 2*nx
   end function se_get_size

   subroutine se_scal(self, alpha)
      class(state_error_vector), intent(inout) :: self
      real(dp),            intent(in)    :: alpha
      call self%state%scal(alpha)
      call self%error%scal(alpha)
   end subroutine se_scal

   subroutine se_axpby(alpha, vec, beta, self)
      class(state_error_vector),        intent(inout) :: self
      class(abstract_vector_rdp), intent(in)    :: vec
      real(dp),                   intent(in)    :: alpha, beta
      select type(vec)
      type is(state_error_vector)
         call self%state%axpby(alpha, vec%state, beta)
         call self%error%axpby(alpha, vec%error, beta)
      class default
         call type_error('vec', 'state_error_vector', 'IN', this_module, 'axpby')
      end select
   end subroutine se_axpby

   subroutine se_rand(self, ifnorm)
      class(state_error_vector), intent(inout) :: self
      logical, optional,   intent(in)    :: ifnorm
      ! internals
      logical :: normalize
      real(dp) :: alpha
      normalize = optval(ifnorm,.true.)
      call self%state%rand(ifnorm=.false.)
      call self%error%rand(ifnorm=.false.)
      if (normalize) then
         alpha = self%state%norm() + self%error%norm()
         call self%state%scal(1.0_dp/alpha)
         call self%error%scal(1.0_dp/alpha)
      endif
   end subroutine se_rand

   !------------------------------------------------------
   !-----     TYPE BOUND PROCEDURES FOR LR STATES    -----
   !------------------------------------------------------

   subroutine initialize_LR_state(self, U, S, rk, rkmax, if_rank_adaptive)
      class(LR_state),            intent(inout) :: self
      type(state_vector),         intent(in)    :: U(:)
      real(dp),                   intent(in)    :: S(:,:)
      integer,                    intent(in)    :: rk
      integer, optional,          intent(in)    :: rkmax
      logical, optional,          intent(in)    :: if_rank_adaptive
      logical                                   :: ifrk

      ! internals
      character(len=*), parameter :: this_procedure = 'initialize_LR_state'
      class(abstract_vector_rdp), allocatable   :: Utmp(:)
      integer :: i, m, rka, n_rem, m_init, info
      character(len=128) :: msg

      ifrk = optval(if_rank_adaptive, .false.)

      ! set time
      call self%reset(full=.true.)

      m = size(U)
      call assert_shape(S, [m,m], 'S', this_module, this_procedure)
      ! optional size argument
      if (present(rkmax)) then
         if (rkmax < rk) then
            call stop_error('rkmax < rk!', this_module, this_procedure)
         end if
         self%rk = rk
         rka = rkmax
         if (ifrk) then
            if (rkmax==rk) then
               call stop_error('rkmax must be larger than rk for rank-adaptive DLRA!', this_module, this_procedure)
            end if
            write(msg,'(A,I0,A)') 'Rank-adaptivity enabled. Computation will begin with X%rk = ', self%rk+1, '.'
            call log_information(msg, this_module, this_procedure)
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
      if (allocated(self%U)) deallocate(self%U)
      if (allocated(self%S)) deallocate(self%S)
      allocate(self%U(rka), source=U(1)); call zero_basis(self%U)
      allocate(self%S(rka,rka)); self%S = 0.0_dp
      write(msg,'(3(A,I0),A)') 'size(X%U) = [ ', rka,' ], X%rk = ', self%rk, ', size(U0) = [ ', m,' ]'
      call log_information(msg, this_module, this_procedure)
      ! copy inputs
      if (self%rk > m) then   ! copy the full IC into self%U
         call copy(self%U(:m), U)
         self%S(:m,:m) = S
         write(msg,'(4X,A,I0,A)') 'Transfer all ', m, ' columns of U0 to X%U.'
         call log_information(msg, this_module, this_procedure)
      else  ! fill the first self%rk columns of self%U with the first self%rk columns of the IC
         call copy(self%U(:self%rk), U(:self%rk))
         self%S(:self%rk,:self%rk) = S(:self%rk,:self%rk)
         write(msg,'(4X,A,I0,A)') 'Transfer the first ', self%rk, ' columns of U0 to X%U.'
         call log_information(msg, this_module, this_procedure)
      end if        

      ! top up basis (to rka for rank-adaptivity) with orthonormal columns if needed
      m_init = min(self%rk, m)
      n_rem = rka - m_init
      if (m > 0) then
         write(msg,'(4X,A,I0,A)') 'Fill remaining ', n_rem, ' columns with orthonormal noise orthonormal to X%U.'
         call log_information(msg, this_module, this_procedure)
         allocate(Utmp(n_rem), source=U(1))
         call initialize_random_orthonormal_basis(Utmp)
         call orthogonalize_against_basis(Utmp, self%U(:m_init), info)
         call copy(self%U(m_init+1:), Utmp)
      end if
   end subroutine initialize_LR_state

end module Ginzburg_Landau_Base
