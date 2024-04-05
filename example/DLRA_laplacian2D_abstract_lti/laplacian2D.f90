module Laplacian2D_LTI
   !> RKLIB module for time integration.
   use rklib_module
   !> exmplib module for exponential propagator
   use LightKrylov_expmlib
   !> LightKrylov for linear algebra.
   use LightKrylov
   use LightROM_AbstractLTIsystems
   !> Standard Library.
   use stdlib_math, only : linspace
   use stdlib_optval, only : optval
   implicit none

   private
   ! problem parameters
   public :: N, nx, dx, dx2, L, rk_b, B, BBT
   ! mesh and operator
   public :: initialize_mesh, laplacian

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
      procedure, pass(self), public :: zero
      procedure, pass(self), public :: dot
      procedure, pass(self), public :: scal
      procedure, pass(self), public :: axpby
      procedure, pass(self), public :: rand
   end type state_vector

   type(state_vector)       :: B(rk_b)
   real(kind=wp)            :: BBT(N**2)

   !-----------------------------------
   !-----     LAPLACE OPERATOR    -----
   !-----------------------------------

   type, extends(abstract_linop), public :: laplace_operator
   contains
      private
      procedure, pass(self), public :: matvec  => direct_matvec_laplace
      procedure, pass(self), public :: rmatvec => direct_matvec_laplace     ! dummy since Lyapunov equation for Laplacian is symmetric
   end type laplace_operator

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

   !------------------------------------------
   !-----     EXPONENTIAL PROPAGATOR     -----
   !------------------------------------------

   type, extends(abstract_linop), public :: rklib_exptA_laplacian
      real(kind=wp) :: tau ! Integration time.
   contains
      private
      procedure, pass(self), public :: matvec  => direct_solver_vec
      procedure, pass(self), public :: rmatvec => direct_solver_vec                  ! dummy
   end type rklib_exptA_laplacian

contains

   !-----     TYPE-BOUND PROCEDURE FOR VECTORS     -----

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

   !-----     TYPE-BOUND PROCEDURE FOR LAPLACE OPERATOR    -----

   subroutine direct_matvec_laplace(self, vec_in, vec_out)
      !> Linear Operator.
      class(laplace_operator),intent(in)  :: self
      !> Input vector.
      class(abstract_vector) , intent(in)  :: vec_in
      !> Output vector.
      class(abstract_vector) , intent(out) :: vec_out
      select type(vec_in)
      type is (state_vector)
         select type(vec_out)
         type is (state_vector)
            call laplacian(vec_out%state, vec_in%state)
         end select
      end select
      return
   end subroutine direct_matvec_laplace

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

   !---------------------------
   !-----    Laplacian    -----
   !---------------------------

   subroutine laplacian(vec_out, vec_in)
      
      !> State vector.
      real(kind=wp)  , dimension(:), intent(in)  :: vec_in
      !> Time-derivative.
      real(kind=wp)  , dimension(:), intent(out) :: vec_out

      !> Internal variables.
      integer             :: i, j, in
      
      in = 1
      vec_out(in)       = (                                  - 4*vec_in(in) + vec_in(in + 1) + vec_in(in + nx)) / dx2
      do in = 2, nx - 1
         vec_out(in)    = (                   vec_in(in - 1) - 4*vec_in(in) + vec_in(in + 1) + vec_in(in + nx)) / dx2
      end do
      in = nx
      vec_out(in)       = (                   vec_in(in - 1) - 4*vec_in(in)                  + vec_in(in + nx)) / dx2
      !
      do i = 2, nx-1
         in = (i-1)*nx + 1
         vec_out(in)    = ( vec_in(in - nx)                  - 4*vec_in(in) + vec_in(in + 1) + vec_in(in + nx)) / dx2
         do j = 2, nx - 1
            in = (i-1)*nx + j
            vec_out(in) = ( vec_in(in - nx) + vec_in(in - 1) - 4*vec_in(in) + vec_in(in + 1) + vec_in(in + nx)) / dx2 
         end do
         in = (i-1)*nx + nx
         vec_out(in)    = ( vec_in(in - nx) + vec_in(in - 1) - 4*vec_in(in)                  + vec_in(in + nx)) / dx2
      end do
      !
      in = N - nx + 1
      vec_out(in)       = ( vec_in(in - nx)                  - 4*vec_in(in) + vec_in(in + 1)                  ) / dx2
      do in = N - nx + 2, N - 1
         vec_out(in)    = ( vec_in(in - nx) + vec_in(in - 1) - 4*vec_in(in) + vec_in(in + 1)                  ) / dx2
      end do
      in = N
      vec_out(in)       = ( vec_in(in - nx) + vec_in(in - 1) - 4*vec_in(in)                                   ) / dx2
         
      return
   end subroutine laplacian

   !-------------------------------------
   !-----    Laplacian for RKlib    -----
   !-------------------------------------

   subroutine rhs(me, t, x, f)
      !> Time-integrator.
      class(rk_class), intent(inout)             :: me
      !> Current time.
      real(kind=wp)  , intent(in)                :: t
      !> State vector.
      real(kind=wp)  , dimension(:), intent(in)  :: x
      !> Time-derivative.
      real(kind=wp)  , dimension(:), intent(out) :: f

      f = 0.0_wp
      call laplacian(f(1:N), x(1:N))
      
      return
   end subroutine rhs

   !-----------------------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE EXPONENTIAL PROPAGATOR RKLIB    -----
   !-----------------------------------------------------------------------------

   subroutine direct_solver_vec(self, vec_in, vec_out)
      !> Linear Operator.
      class(rklib_exptA_laplacian), intent(in)  :: self
      !> Input vector.
      class(abstract_vector) , intent(in)  :: vec_in
      !> Output vector.
      class(abstract_vector) , intent(out) :: vec_out

      !> Time-integrator.
      type(rks54_class) :: prop
      real(kind=wp)     :: dt = 1.0_wp

      select type(vec_in)
      type is (state_vector)
         select type(vec_out)
         type is (state_vector)

             !> Initialize propagator.
             call prop%initialize(n=N, f=rhs)
             !> Integrate forward in time.
             call prop%integrate(0.0_wp, vec_in%state, dt, self%tau, vec_out%state)

         end select
      end select
      return
   end subroutine direct_solver_vec

end module Laplacian2D_LTI