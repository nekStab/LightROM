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

end module laplacian2D_LTI_Lyapunov_Base