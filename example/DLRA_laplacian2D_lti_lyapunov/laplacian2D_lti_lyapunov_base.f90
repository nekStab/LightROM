module laplacian2D_LTI_Lyapunov_Base
   !> LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov_utils, only : assert_shape
   use LightROM_AbstractLTIsystems
   !> Standard Library.
   use stdlib_math, only : linspace
   use stdlib_optval, only : optval
   implicit none

   private
   ! problem parameters
   public :: N, nx, dx, dx2, L, rk_b, B, BBT
   ! mesh and operator
   public :: initialize_mesh
   ! utils
   public :: get_state, set_state, init_rand

   !------------------------------
   !-----     PARAMETERS     -----
   !------------------------------

   ! --> Mesh related parameters.
   real(kind=wp), parameter :: L  = 1.0_wp  !> Domain length
   integer,       parameter :: nx = 4      !> Number of grid points per direction
   integer,       parameter :: N  = nx**2   !> total number of grid points
   real(kind=wp), parameter :: dx = L/nx    !> Grid size.
   real(kind=wp), parameter :: dx2= dx**2   !> Grid size.
   integer,       parameter :: rk_b = 5     !> rank of the RHS

   !-------------------------------------------
   !-----     LIGHTKRYLOV VECTOR TYPE     -----
   !-------------------------------------------

   type, extends(abstract_vector), public :: state_vector
      real(kind=wp) :: state(N) = 0.0_wp
      contains
      private
      procedure, pass(self), public :: zero => vector_zero
      procedure, pass(self), public :: dot => vector_dot
      procedure, pass(self), public :: scal => vector_scal
      procedure, pass(self), public :: axpby => vector_axpby
      procedure, pass(self), public :: rand => vector_rand
   end type state_vector

   !-------------------------------------------
   !-----     LIGHTKRYLOV VECTOR TYPE     -----
   !-------------------------------------------

   type, extends(abstract_vector), public :: state_matrix
      real(kind=wp) :: state(N**2) = 0.0_wp
   contains
      private
      procedure, pass(self), public :: zero => matrix_zero
      procedure, pass(self), public :: dot => matrix_dot
      procedure, pass(self), public :: scal => matrix_scal
      procedure, pass(self), public :: axpby => matrix_axpby
      procedure, pass(self), public :: rand => matrix_rand
   end type state_matrix

   type(state_vector)       :: B(rk_b)
   real(kind=wp)            :: BBT(N**2)

   !-----------------------------------------------
   !-----     LIGHTKRYLOV LTI SYSTEM TYPE     -----
   !-----------------------------------------------

   type, extends(abstract_lti_system), public :: lti_system
   end type lti_system

   !-------------------------------------------------------
   !-----     LIGHTKRYLOV SYM LOW RANK STATE TYPE     -----
   !-------------------------------------------------------

   type, extends(abstract_sym_low_rank_state), public :: LR_state
   contains
      private
      procedure, pass(self), public :: set_LR_state
   end type LR_state

contains

   !-----     TYPE-BOUND PROCEDURE FOR VECTORS     -----

   subroutine vector_zero(self)
      class(state_vector), intent(inout) :: self
      self%state = 0.0_wp
      return
   end subroutine vector_zero

   real(kind=wp) function vector_dot(self, vec) result(alpha)
      class(state_vector)   , intent(in) :: self
      class(abstract_vector), intent(in) :: vec
      select type(vec)
      type is (state_vector)
         alpha = dot_product(self%state, vec%state)
      end select
      return
   end function vector_dot

   subroutine vector_scal(self, alpha)
      class(state_vector), intent(inout) :: self
      real(kind=wp)      , intent(in)    :: alpha
      self%state = self%state * alpha
      return
   end subroutine vector_scal

   subroutine vector_axpby(self, alpha, vec, beta)
      class(state_vector)   , intent(inout) :: self
      class(abstract_vector), intent(in)    :: vec
      real(kind=wp)         , intent(in)    :: alpha, beta
      select type(vec)
      type is (state_vector)
         self%state = alpha*self%state + beta*vec%state
      end select
      return
   end subroutine vector_axpby

   subroutine vector_rand(self, ifnorm)
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
   end subroutine vector_rand

   !-----     TYPE-BOUND PROCEDURE FOR MATRICES     -----

   subroutine matrix_zero(self)
      class(state_matrix), intent(inout) :: self
      self%state = 0.0_wp
      return
   end subroutine matrix_zero

   real(kind=wp) function matrix_dot(self, vec) result(alpha)
      class(state_matrix)   , intent(in) :: self
      class(abstract_vector), intent(in) :: vec
      select type(vec)
      type is(state_matrix)
          alpha = dot_product(self%state, vec%state)
      end select
      return
   end function matrix_dot

   subroutine matrix_scal(self, alpha)
      class(state_matrix), intent(inout) :: self
      real(kind=wp)      , intent(in)    :: alpha
      self%state = self%state * alpha
      return
   end subroutine matrix_scal  

   subroutine matrix_axpby(self, alpha, vec, beta)
      class(state_matrix)   , intent(inout) :: self
      class(abstract_vector), intent(in)    :: vec
      real(kind=wp)         , intent(in)    :: alpha, beta
      select type(vec)
      type is(state_matrix)
          self%state = alpha*self%state + beta*vec%state
      end select
      return
   end subroutine matrix_axpby

   subroutine matrix_rand(self, ifnorm)
      class(state_matrix), intent(inout) :: self
      logical, optional,   intent(in)    :: ifnorm
      ! internals
      logical :: normalize
      real(kind=wp) :: alpha
      normalize = optval(ifnorm, .true.)
      call random_number(self%state)
      if (normalize) then
         alpha = self%norm()
         call self%scal(1.0/alpha)
      endif
      return
   end subroutine matrix_rand

   !---------------------------------------
   !-----     CONSTRUCT THE MESH      -----
   !---------------------------------------

   subroutine initialize_mesh()
      implicit none
      !> Mesh array.
      real(kind=wp), allocatable :: x(:)
      integer :: i

      !> Construct mesh.
      x = linspace(-L/2, L/2, nx)

      return
   end subroutine initialize_mesh

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
      type is (state_matrix)
         call assert_shape(mat_out, (/ N, N /), 'get_state -> state_matrix', 'mat_out')
         mat_out = reshape(state_in(1)%state, (/ N, N /))
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
      type is (state_matrix)
         call assert_shape(mat_in, (/ N, N /), 'set_state -> state_matrix', 'mat_in')
         call mat_zero(state_out)
         state_out(1)%state = reshape(mat_in, shape(state_out(1)%state))
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
      type is (state_matrix)
         kdim = size(state)
         do k = 1, kdim
            call state(k)%rand(ifnorm = normalize)
         end do
      end select
      return
   end subroutine init_rand

   !------------------------------------------------------------
   !-----     UTILITIES FOR SYM LOW RANK REPRESENTATION    -----
   !------------------------------------------------------------

   subroutine set_LR_state(self, U, S)
      class(LR_state), intent(inout) :: self
      real(kind=wp),   intent(in)    :: U(:,:)
      real(kind=wp),   intent(in)    :: S(:,:)
      ! internals
      integer :: rk
      rk = size(U,2)
      call assert_shape(S, (/ rk, rk /), 'set_LR_state', 'S')
      call set_state(self%U, U)
      self%S = S
      return
   end subroutine set_LR_state

end module laplacian2D_LTI_Lyapunov_Base