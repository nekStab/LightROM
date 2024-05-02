module Ginzburg_Landau_Base
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov_utils, only : assert_shape
   ! LightROM
   use LightROM_AbstractLTIsystems
   ! Standard Library.
   use stdlib_math, only : linspace
   use stdlib_optval, only : optval
   implicit none

   private
   public :: nx, dx
   public :: nu, gamma, mu_0, c_mu, mu_2, mu
   public :: rk_b, x_b, s_b, rk_c, x_c, s_c
   public :: B, CT, weight
   public :: initialize_parameters
   public :: set_state, get_state, init_rand

   !-------------------------------
   !-----     PARAMETERS 1    -----
   !-------------------------------

   ! Mesh related parameters.
   real(kind=wp), parameter :: L  = 50.0_wp ! Domain length
   integer      , parameter :: nx = 128     ! Number of grid points (excluding boundaries).
   real(kind=wp)            :: dx           ! Grid size.

   !-------------------------------------------
   !-----     LIGHTKRYLOV VECTOR TYPE     -----
   !-------------------------------------------

   type, extends(abstract_vector), public :: state_vector
      real(kind=wp) :: state(2*nx) = 0.0_wp
   contains
      private
      procedure, pass(self), public :: zero
      procedure, pass(self), public :: dot
      procedure, pass(self), public :: scal
      procedure, pass(self), public :: axpby
      procedure, pass(self), public :: rand
   end type state_vector

   !-------------------------------------------------------
   !-----     LIGHTKRYLOV SYM LOW RANK STATE TYPE     -----
   !-------------------------------------------------------

   type, extends(abstract_sym_low_rank_state), public :: LR_state
   contains
      private
      procedure, pass(self), public :: initialize_LR_state
   end type LR_state

   !-------------------------------
   !-----     PARAMETERS 2    -----
   !-------------------------------

   ! Physical parameters.
   complex(kind=wp), parameter :: nu    = cmplx(2.0_wp, 0.2_wp, kind=wp)
   complex(kind=wp), parameter :: gamma = cmplx(1.0_wp, -1.0_wp, kind=wp)
   real(kind=wp)   , parameter :: mu_0  = 0.38_wp
   real(kind=wp)   , parameter :: c_mu  = 0.2_wp
   real(kind=wp)   , parameter :: mu_2  = -0.01_wp
   real(kind=wp)               :: mu(1:nx)

   ! Input-Output system parameters
   real(kind=wp)               :: weight(2*nx)       ! integration weights
   integer,       parameter    :: rk_b = 1           ! number of inputs to the system
   real(kind=wp), parameter    :: x_b = -11.0_wp     ! location of input Gaussian
   real(kind=wp), parameter    :: s_b = 1.0_wp       ! variance of input Gaussian
   type(state_vector)          :: B(rk_b)
   real(kind=wp), parameter    :: x_c = sqrt(-2.0_wp*(mu_0 - c_mu**2)/mu_2) ! location of input Gaussian
   real(kind=wp), parameter    :: s_c = 1.0_wp       ! variance of input Gaussian
   integer,       parameter    :: rk_c = 1           ! number of outputs to the system
   type(state_vector)          :: CT(rk_c)

contains

   !--------------------------------------------------------------
   !-----     CONSTRUCT THE MESH AND PHYSICAL PARAMETERS     -----
   !--------------------------------------------------------------

   subroutine initialize_parameters()
      implicit none
      ! Mesh array.
      real(kind=wp), allocatable :: x(:)
      real(kind=wp)              :: x2(1:2*nx)

      ! Construct mesh.
      x = linspace(-L/2, L/2, nx+2)
      dx = x(2)-x(1)

      ! Construct mu(x)
      mu(:) = (mu_0 - c_mu**2) + (mu_2 / 2.0_wp) * x(2:nx+1)**2

      ! Define integration weights
      weight       = dx
      !weight(1)    = 0.5_wp*dx
      !weight(nx)   = 0.5_wp*dx
      !weight(nx+1) = 0.5_wp*dx
      !weight(2*nx) = 0.5_wp*dx

      ! Construct B & C
      x2(1:nx)      = x(2:nx+1)
      x2(nx+1:2*nx) = x(2:nx+1)
      ! actuator is a Guassian centered just upstream of branch I
      B(1)%state = exp(-((x2 - x_b)/s_b)**2)
      ! the sensor is a Gaussian centered at branch II
      CT(1)%state = exp(-((x2 - x_c)/s_c)**2)      

      return
   end subroutine initialize_parameters

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

   real(kind=wp) function dot(self, vec) result(alpha)
      ! weighted inner product
      class(state_vector)   , intent(in) :: self
      class(abstract_vector), intent(in) :: vec
      select type(vec)
      type is(state_vector)
         alpha = dot_product(self%state, weight*vec%state)
      end select
      return
   end function dot

   subroutine scal(self, alpha)
      class(state_vector), intent(inout) :: self
      real(kind=wp)      , intent(in)    :: alpha
      self%state = self%state * alpha
      return
   end subroutine scal

   subroutine axpby(self, alpha, vec, beta)
      class(state_vector)   , intent(inout) :: self
      class(abstract_vector), intent(in)    :: vec
      real(kind=wp)         , intent(in)    :: alpha, beta
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
      real(kind=wp) :: alpha
      normalize = optval(ifnorm,.true.)
      call random_number(self%state)
      if (normalize) then
         alpha = self%norm()
         call self%scal(1.0/alpha)
      endif
      return
   end subroutine rand

   !----------------------------------------------
   !-----     UTILITIES FOR STATE_VECTORS    -----
   !----------------------------------------------

   subroutine get_state(mat_out, state_in)
      !! Utility function to transfer data from a state vector to a real array
      real(kind=wp),          intent(out) :: mat_out(:,:)
      class(abstract_vector), intent(in)  :: state_in(:)
      ! internal variables
      integer :: k, kdim
      mat_out = 0.0_wp
      select type (state_in)
      type is (state_vector)
         kdim = size(state_in)
         call assert_shape(mat_out, (/ 2*nx, kdim /), 'get_state -> state_vector', 'mat_out')
         do k = 1, kdim
            mat_out(:,k) = state_in(k)%state
         end do
      end select
      return
   end subroutine get_state

   subroutine set_state(state_out, mat_in)
      !! Utility function to transfer data from a real array to a state vector
      class(abstract_vector), intent(out) :: state_out(:)
      real(kind=wp),          intent(in)  :: mat_in(:,:)
      ! internal variables
      integer       :: k, kdim
      select type (state_out)
      type is (state_vector)
         kdim = size(state_out)
         call assert_shape(mat_in, (/ 2*nx, kdim /), 'set_state -> state_vector', 'mat_in')
         call mat_zero(state_out)
         do k = 1, kdim
            state_out(k)%state = mat_in(:,k)
         end do
      end select
      return
   end subroutine set_state

   subroutine init_rand(state, ifnorm)
      !! Utility function to initialize a state vector with random data
      class(abstract_vector), intent(inout)  :: state(:)
      logical, optional,      intent(in)     :: ifnorm
      ! internal variables
      integer :: k, kdim
      logical :: normalize
      normalize = optval(ifnorm,.true.)
      select type (state)
      type is (state_vector)
         kdim = size(state)
         do k = 1, kdim
            call state(k)%rand(ifnorm = normalize)
         end do
      end select
      return
   end subroutine init_rand

   !------------------------------------------------------
   !-----     TYPE BOUND PROCEDURES FOR LR STATES    -----
   !------------------------------------------------------

   subroutine initialize_LR_state(self, U, S, rk)
      class(LR_state),        intent(inout) :: self
      class(abstract_vector), intent(in)    :: U(:)
      real(kind=wp),          intent(in)    :: S(:,:)
      integer,                intent(in)    :: rk

      if (rk > size(U)) then
         write(*,*) 'Input state rank is lower than the chosen rank! Abort.'
         STOP 1
         ! this could be improved by initialising extra columns with random vectors
         ! orthonormalize these against the existing columns of U and set the corresponding
         ! entries in S to 0.
      end if

      select type (U)
      type is (state_vector)
         allocate(self%U(1:rk), source=U(1:rk))
         allocate(self%S(1:rk,1:rk)); 
         self%S(1:rk,1:rk) = S(1:rk,1:rk) 
      end select
      return
   end subroutine initialize_LR_state

end module Ginzburg_Landau_Base
