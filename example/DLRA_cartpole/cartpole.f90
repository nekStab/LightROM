module cartpole
   use LightKrylov
   use lightKrylov_utils, only: assert_shape
   use rklib_module
   use stdlib_linalg, only: eye
   use stdlib_optval, only: optval

   implicit none

   private
   ! problem parameters
   public :: n, m, p, delta, g, m_pole, m_cart, L, A, B, Ctrans, BRinvBT, CTQC, Q, R
   ! setup/problem subroutines
   public :: initialize_cartpole

   integer,       parameter ::      n = 4        ! number of degrees of freedom
   integer,       parameter ::      m = 1        ! number of actuators
   integer,       parameter ::      p = 4        ! number of observables
   real(kind=wp), parameter ::  delta =   1.0_wp ! friction coefficient
   real(kind=wp), parameter ::      g = -10.0_wp ! gravitational acceleration
   real(kind=wp), parameter :: m_pole =   1.0_wp ! mass of the pole
   real(kind=wp), parameter :: m_cart =   5.0_wp ! mass of the cart
   real(kind=wp), parameter ::      L =   2.0_wp ! length of the pole
   real(kind=wp)            :: A(n,n)       ! system matrix
   real(kind=wp)            :: R(m,m)       ! control weights
   real(kind=wp)            :: Q(n,n)       ! state weights
   real(kind=wp)            :: BRinvBT(n,n) ! -P
   real(kind=wp)            :: CTQC(n,n)    ! Q

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
   type(state_vector) :: Ctrans(p)

   type, extends(abstract_linop), public :: cartpole_op
   contains
      private
      procedure, pass(self), public :: matvec  => direct_matvec_cp
      procedure, pass(self), public :: rmatvec => direct_matvec_cp ! dummy (not used)
   end type cartpole_op

   interface apply_Q
      module procedure apply_Q_mat
      module procedure apply_Q_state
   end interface apply_Q

   interface apply_Rinv
      module procedure apply_Rinv_mat
      module procedure apply_Rinv_state
   end interface apply_Rinv

contains

   subroutine initialize_cartpole
      implicit none

      real(kind=wp) :: Id(n,n)
      integer       :: i
      ! reference
      real(kind=wp) :: Bdata(n,m)
      real(kind=wp) :: Bwrk(n,m)
      real(kind=wp) :: Cdata(p,n)
      real(kind=wp) :: Cwrk(p,n)

      Id = eye(n)

      A = transpose(reshape((/ 0.0_wp,            1.0_wp,                          0.0_wp, 0.0_wp, &
                             & 0.0_wp,     -delta/m_cart,                 m_pole*g/m_cart, 0.0_wp, &
                             & 0.0_wp,            0.0_wp,                          0.0_wp, 1.0_wp, &
                             & 0.0_wp, -delta/(m_cart*L), -(m_pole + m_cart)*g/(m_cart*L), 0.0_wp /), shape(A)))

      B(1)%state = (/ 0.0_wp, 1.0/m_cart, 0.0_wp, 1.0/(m_cart*L) /)

      do i = 1, p
         Ctrans(i)%state = Id(:,i)
      end do

      Bdata = 0.0_wp; Bwrk  = 0.0_wp; Cdata = 0.0_wp; Cwrk = 0.0_wp
      call get_state(Bdata, B);      call get_state(Bwrk,  B)
      call get_state(Cdata, Ctrans); call get_state(Cwrk,  Ctrans)
      call apply_Rinv(Bwrk)
      BRinvBT = matmul( Bwrk, transpose(Bdata) )
      call apply_Q(Cwrk)
      CTQC    = matmul( Cwrk, transpose(Cdata) )
      
      return

   end subroutine initialize_cartpole

   subroutine apply_Q_mat(X)
      implicit none
      real(kind=wp), intent(inout) :: X(:,:)
      !
      ! Q = eye(n) --> could also just return immediately
      !
      ! internals
      real(kind=wp), parameter :: Qdiag(1:p) = (/ 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp /)
      integer :: i

      if (size(X,2) /= p) then
         write (*,*) "INFO: In apply_Q_mat input matrix X has the wrong size"
         write (*,*) "Expected shape is ", n,"x",p
         call stop_error("Aborting due to illegal matrix operation.")
         STOP 1
      end if
      
      do i = 1, size(X,2)
         X(:,i) = X(:,i)*Qdiag(i)
      end do
      return
   end subroutine apply_Q_mat

   subroutine apply_Q_state(X)
      implicit none
      class(abstract_vector), intent(inout) :: X(:)
      !
      ! Q = eye(n) --> could also just return immediately
      !
      ! internals
      real(kind=wp), parameter :: Qdiag(1:p) = (/ 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp /)
      integer :: i

      if (size(X) /= p) then
         write (*,*) "INFO: In apply_Q_state input basis X has the wrong size"
         write (*,*) "Expected shape is ",p
         call stop_error("Aborting due to illegal matrix operation.")
         STOP 1
      end if

      select type (X)
      type is (state_vector)
         do i = 1, size(X)
            call X(i)%scal(Qdiag(i)) 
         end do
      end select
      return
   end subroutine apply_Q_state

   subroutine apply_Rinv_mat(U)
      implicit none
      real(kind=wp), intent(inout) :: U(:,:)
      !
      ! R = eye(m) --> could also just return immediately
      !
      ! internals
      real(kind=wp), parameter :: Rdiag(1) = (/ 0.001_wp /)
      integer :: i
      do i = 1, size(U,2)
         U(:,i) = U(:,i)/Rdiag(i) 
      end do
      return
   end subroutine apply_Rinv_mat

   subroutine apply_Rinv_state(U)
      implicit none
      class(abstract_vector), intent(inout) :: U(:)
      !
      ! R = eye(m) --> could also just return immediately
      !
      ! internals
      real(kind=wp), parameter :: Rdiag(1) = (/ 1.0_wp /)
      integer :: i
      select type (U)
      type is (state_vector)
         do i = 1, size(U)
            call U(i)%scal(1.0/Rdiag(i)) 
         end do
      end select
      return
   end subroutine apply_Rinv_state

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
      class(cartpole_op),intent(in)        :: self
      !> Input vector.
      class(abstract_vector) , intent(in)  :: vec_in
      !> Output vector.
      class(abstract_vector) , intent(out) :: vec_out

      select type(vec_in)
      type is (state_vector)
         select type(vec_out)
         type is (state_vector)
            vec_out%state = matmul(A, vec_in%state)
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