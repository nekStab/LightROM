module cartpole
   use LightKrylov
   use lightrom_AbstractLTIsystems
   use lightKrylov_utils, only: assert_shape
   use rklib_module
   use stdlib_linalg, only: eye
   use stdlib_optval, only: optval

   implicit none

   private
   ! problem parameters
   public :: n, m, delta, g, m_pole, m_cart, L, Amat, B, BRinvBTmat, Qmat, Rinv
   ! setup/problem subroutines
   public :: initialize_cartpole, get_state, set_state

   integer,       parameter ::      n = 4        ! number of degrees of freedom
   integer,       parameter ::      m = 1        ! number of actuators
   real(kind=wp), parameter ::  delta =   1.0_wp ! friction coefficient
   real(kind=wp), parameter ::      g = -10.0_wp ! gravitational acceleration
   real(kind=wp), parameter :: m_pole =   1.0_wp ! mass of the pole
   real(kind=wp), parameter :: m_cart =   5.0_wp ! mass of the cart
   real(kind=wp), parameter ::      L =   2.0_wp ! length of the pole
   real(kind=wp)            :: Amat(n,n)         ! system matrix
   real(kind=wp)            :: BRinvBTmat(n,n)   ! -P
   real(kind=wp)            :: Qmat(n,n)         ! Q
   real(kind=wp)            :: Rinv(m,m)

   !-------------------------------------------
   !-----     LIGHTKRYLOV VECTOR TYPE     -----
   !-------------------------------------------

   type, extends(abstract_vector), public :: state_vector
      real(kind=wp) :: state(n) = 0.0_wp         ! x, \dot{x}, \theta, \dot{\theta}
      contains
      procedure, pass(self), public :: zero
      procedure, pass(self), public :: dot
      procedure, pass(self), public :: scal
      procedure, pass(self), public :: axpby
      procedure, pass(self), public :: rand
   end type state_vector

   type(state_vector) :: B(m)

   !-----------------------------------------------
   !-----     LIGHTKRYLOV LTI SYSTEM TYPE     -----
   !-----------------------------------------------

   type, extends(abstract_lti_system), public :: lti_system
   end type lti_system

   type, extends(abstract_linop), public :: cartpole_linop
   contains
      private
      procedure, pass(self), public :: matvec  => direct_matvec_cp
      procedure, pass(self), public :: rmatvec => direct_matvec_cp ! dummy (not used)
   end type cartpole_linop

contains

   subroutine initialize_cartpole
      implicit none

      ! reference
      real(kind=wp) :: Bdata(n,m)
      real(kind=wp) :: Bwrk(n,m)

      Amat = transpose(reshape((/ 0.0_wp,            1.0_wp,                          0.0_wp, 0.0_wp, &
                                & 0.0_wp,     -delta/m_cart,                 m_pole*g/m_cart, 0.0_wp, &
                                & 0.0_wp,            0.0_wp,                          0.0_wp, 1.0_wp, &
                                & 0.0_wp, -delta/(m_cart*L), -(m_pole + m_cart)*g/(m_cart*L), 0.0_wp /), shape(Amat)))

      B(1)%state = (/ 0.0_wp, 1.0/m_cart, 0.0_wp, 1.0/(m_cart*L) /)

      Rinv = 1e3*eye(m)

      Bdata = 0.0_wp; Bwrk  = 0.0_wp;
      call get_state(Bwrk, B);  
      Bdata = matmul(Bwrk, Rinv)
      call get_state(Bwrk, B)
      BRinvBTmat = matmul( Bdata, transpose(Bwrk) )
      
      Qmat = eye(n)
      
      return

   end subroutine initialize_cartpole

   !----------------------------------------------------
   !-----     TYPE-BOUND PROCEDURE FOR VECTORS     -----
   !----------------------------------------------------

   subroutine zero(self)
      class(state_vector), intent(inout) :: self
      self%state = 0.0_wp
      return
   end subroutine zero

   real(kind=wp) function dot(self, vec) result(alpha)
      class(state_vector)   , intent(in) :: self
      class(abstract_vector), intent(in) :: vec
      select type(vec)
      type is (state_vector)
         alpha = dot_product(self%state, vec%state)
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
      type is (state_vector)
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

   !--------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR CARTPOLE OPERATOR    -----
   !--------------------------------------------------------------

   subroutine direct_matvec_cp(self, vec_in, vec_out)
      !> Linear Operator.
      class(cartpole_linop),intent(in)     :: self
      !> Input vector.
      class(abstract_vector) , intent(in)  :: vec_in
      !> Output vector.
      class(abstract_vector) , intent(out) :: vec_out

      select type(vec_in)
      type is (state_vector)
         select type(vec_out)
         type is (state_vector)
            vec_out%state = matmul(Amat, vec_in%state)
         end select
      end select
   
      return
   end subroutine direct_matvec_cp

   !--------------------------------------------------------------------
   !-----     UTILITIES FOR STATE_VECTOR AND STATE MATRIX TYPES    -----
   !--------------------------------------------------------------------

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
         call assert_shape(mat_out, (/ N, kdim /), 'get_state -> state_vector', 'mat_out')
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
         call assert_shape(mat_in, (/ N, kdim /), 'set_state -> state_vector', 'mat_in')
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

end module cartpole