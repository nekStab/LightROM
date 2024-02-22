module Laplacian
  !> RKLIB module for time integration.
  use rklib_module
  !> LightKrylov for linear algebra.
  use LightKrylov
  !use LightROM
  !> Standard Library.
  use stdlib_math, only : linspace
  implicit none

  private
  public :: nx, dx, L
  public :: initialize_mesh

  !------------------------------
  !-----     PARAMETERS     -----
  !------------------------------

  ! --> Mesh related parameters.
  real(kind=wp), parameter :: L  = 1.0_wp  !> Domain length
  integer      , parameter :: nx = 10      !> Number of grid points
  real(kind=wp), parameter :: dx = L/nx    !> Grid size.

  !-------------------------------------------
  !-----     LIGHTKRYLOV VECTOR TYPE     -----
  !-------------------------------------------

  type, extends(abstract_vector), public :: state_vector
     real(kind=wp) :: state(nx) = 0.0_wp
   contains
     private
     procedure, pass(self), public :: zero
     procedure, pass(self), public :: dot
     procedure, pass(self), public :: scal
     procedure, pass(self), public :: axpby
  end type state_vector

  !------------------------------------------
  !-----     EXPONENTIAL PROPAGATOR     -----
  !------------------------------------------

  type, extends(abstract_linop), public :: exponential_prop
     real(kind=wp) :: tau ! Integration time.
   contains
     private
     procedure, pass(self), public :: matvec => direct_solver
     procedure, pass(self), public :: rmatvec => adjoint_solver
  end type exponential_prop

contains

  !=============================
  !=============================
  !=====                   =====
  !=====     LAPLACIAN     =====
  !=====                   =====
  !=============================
  !=============================

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


  !------------------------------
  !-----      Laplacian     -----
  !------------------------------

  subroutine rhs(me, t, x, f)
    !> Time-integrator.
    class(rk_class), intent(inout)             :: me
    !> Current time.
    real(kind=wp)  , intent(in)                :: t
    !> State vector.
    real(kind=wp)  , dimension(:), intent(in)  :: x
    !> Time-derivative.
    real(kind=wp)  , dimension(:), intent(out) :: f

    !> Internal variables.
    integer :: i, j, k
    real(kind=wp), dimension(nx) :: v, dv

    !> Sets the internal variables.
    v   = x(1:nx)
    f   = 0.0_wp
    dv  = 0.0_wp

    !> Left most boundary points 
    dv(1) = (- 2*v(1) + v(2)) / dx**2   !> Diffusion term (zero Dirichlet)
    !> Interior nodes.
    do i = 2, nx-1
       dv(i) = (v(i-1) - 2*v(i) + v(i+1)) / dx**2
    enddo
    !> Right most boundary points
    dv(nx) = (v(nx-1) - 2*v(nx)) / dx**2

    !> Combine the parts and copy result to the output array.
    f(1:nx) = dv
    
    return
  end subroutine rhs

  !-------------------------------------
  !-----     Adjoint diffusion     -----
  !-------------------------------------

  subroutine adjoint_rhs(me, t, x, f)
    !> Time-integrator.
    class(rk_class), intent(inout)             :: me
    !> Current time.
    real(kind=wp)  , intent(in)                :: t
    !> State vector.
    real(kind=wp)  , dimension(:), intent(in)  :: x
    !> Time-derivative.
    real(kind=wp)  , dimension(:), intent(out) :: f

    ! this is a dummy function since the laplacian is self-adjoint
    f  = 0.0_wp
    call rhs(me, t, x, f)

    return
  end subroutine adjoint_rhs

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
    class(state_vector)   , intent(in) :: self
    class(abstract_vector), intent(in) :: vec
    select type(vec)
    type is(state_vector)
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
    type is(state_vector)
       self%state = alpha*self%state + beta*vec%state
    end select
    return
  end subroutine axpby

  !------------------------------------------------------------------------
  !-----     TYPE-BOUND PROCEDURES FOR THE EXPONENTIAL PROPAGATOR     -----
  !------------------------------------------------------------------------

  subroutine direct_solver(self, vec_in, vec_out)
    !> Linear Operator.
    class(exponential_prop), intent(in)  :: self
    !> Input vector.
    class(abstract_vector) , intent(in)  :: vec_in
    !> Output vector.
    class(abstract_vector) , intent(out) :: vec_out

    !> Time-integrator.
    type(rks54_class) :: prop
    real(kind=wp)     :: dt = 1.0_wp

    select type(vec_in)
    type is(state_vector)
       select type(vec_out)
       type is(state_vector)

          !> Initialize propagator.
          call prop%initialize(n=nx, f=rhs)
          !> Integrate forward in time.
          call prop%integrate(0.0_wp, vec_in%state, dt, self%tau, vec_out%state)

       end select
    end select
    return
  end subroutine direct_solver

  subroutine adjoint_solver(self, vec_in, vec_out)
    !> Linear Operator.
    class(exponential_prop), intent(in)  :: self
    !> Input vector.
    class(abstract_vector) , intent(in)  :: vec_in
    !> Output vector.
    class(abstract_vector) , intent(out) :: vec_out

    !> Time-integrator.
    type(rks54_class) :: prop
    real(kind=wp)     :: dt = 1.0_wp

    select type(vec_in)
    type is(state_vector)
       select type(vec_out)
       type is(state_vector)

          !> Initialize propagator.
          call prop%initialize(n=nx, f=adjoint_rhs)
          !> Integrate forward in time.
          call prop%integrate(0.0_wp, vec_in%state, dt, self%tau, vec_out%state)

       end select
    end select
    return
  end subroutine adjoint_solver

end module Laplacian
