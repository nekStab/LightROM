module Laplacian_Lyapunov
  use Laplacian
  !> RKLIB module for time integration.
  use rklib_module
  !> LightKrylov for linear algebra.
  use LightKrylov
  !use LightROM
  !> Standard Library.
  use stdlib_math, only : linspace
  implicit none

  private
  public :: initialize_mesh_lyap
  
  !-------------------------------------------
  !-----     LIGHTKRYLOV VECTOR TYPE     -----
  !-------------------------------------------

  type, extends(abstract_vector), public :: state_vector_lyap
     real(kind=wp) :: state(nx**2) = 0.0_wp
   contains
     private
     procedure, pass(self), public :: zero
     procedure, pass(self), public :: dot
     procedure, pass(self), public :: scal
     procedure, pass(self), public :: axpby
  end type state_vector_lyap

  !------------------------------------------
  !-----     EXPONENTIAL PROPAGATOR     -----
  !------------------------------------------

  type, extends(abstract_linop), public :: exponential_prop_lyap
     real(kind=wp) :: tau ! Integration time.
   contains
     private
     procedure, pass(self), public :: matvec => direct_solver_lyap
     procedure, pass(self), public :: rmatvec => adjoint_solver_lyap
  end type exponential_prop_lyap

contains

  !=====================================
  !=====================================
  !=====                           =====
  !=====     LAPLACIAN LYAPUNOV    =====
  !=====                           =====
  !=====================================
  !=====================================

  !---------------------------------------
  !-----     CONSTRUCT THE MESH      -----
  !---------------------------------------

  subroutine initialize_mesh_lyap()
    implicit none
    !> Mesh array.
    real(kind=wp), allocatable :: x0(:)
    real(kind=wp), allocatable :: x(:)
    integer :: i

    allocate(x(1:nx**2))

    !> Construct mesh.
    x0 = linspace(-L/2, L/2, nx)
    do i = 1, nx
      x((i-1)*nx+1:i*nx) = x0
    enddo
    
    return
  end subroutine initialize_mesh_lyap


  !------------------------------
  !-----      Laplacian     -----
  !------------------------------

  subroutine rhs_lyap(me, t, x, f)
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
    real(kind=wp), dimension(nx,nx) :: v, dv, dvT
    real(kind=wp), dimension(nx**2) :: q

    !> Sets the internal variables.
    v   = reshape(x(1:nx**2),(/nx, nx/))
    f   = 0.0_wp
    dv  = 0.0_wp
    dvT = 0.0_wp
    q   = 1.0_wp

    !> We compute the action of the Lyapunov operator without using the transpose of A
    !           L(X) = A @ X +   X @ A.T     + Q
    !                = A @ X + ( A @ X.T ).T + Q
    do j = 1,nx
      !> Left most boundary points 
      dv( 1,j) = (- 2*v(1,j) + v(2,j)) / dx**2   !> Diffusion term (zero Dirichlet)
      dvT(j,1) = (- 2*v(j,1) + v(j,2)) / dx**2   !> Diffusion term (zero Dirichlet)
      !> Interior nodes.
      do i = 2, nx-1
         dv( i,j) = (v(i+1,j) - 2*v(i,j) + v(i-1,j)) / dx**2
         dvT(j,i) = (v(j,i+1) - 2*v(j,i) + v(j,i-1)) / dx**2
      enddo
      !> Right most boundary points
      dv( nx,j) = (v(nx-1,j) - 2*v(nx,j)) / dx**2
      dvT(j,nx) = (v(j,nx-1) - 2*v(j,nx)) / dx**2
    enddo

    !> Combine the parts and copy result to the output array.
    f(1:nx**2) = reshape(dv, shape(q)) + reshape(dvT, shape(q)) + q
    
    return
  end subroutine rhs_lyap

  !-------------------------------------
  !-----     Adjoint diffusion     -----
  !-------------------------------------

  subroutine adjoint_rhs_lyap(me, t, x, f)
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
    call rhs_lyap(me, t, x, f)

    return
  end subroutine adjoint_rhs_lyap

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
    class(state_vector_lyap), intent(inout) :: self
    self%state = 0.0_wp
    return
  end subroutine zero

  real(kind=wp) function dot(self, vec) result(alpha)
    class(state_vector_lyap)   , intent(in) :: self
    class(abstract_vector), intent(in) :: vec
    select type(vec)
    type is(state_vector_lyap)
       alpha = dot_product(self%state, vec%state)
    end select
    return
  end function dot

  subroutine scal(self, alpha)
    class(state_vector_lyap), intent(inout) :: self
    real(kind=wp)      , intent(in)    :: alpha
    self%state = self%state * alpha
    return
  end subroutine scal

  subroutine axpby(self, alpha, vec, beta)
    class(state_vector_lyap)   , intent(inout) :: self
    class(abstract_vector), intent(in)    :: vec
    real(kind=wp)         , intent(in)    :: alpha, beta
    select type(vec)
    type is(state_vector_lyap)
       self%state = alpha*self%state + beta*vec%state
    end select
    return
  end subroutine axpby

  !------------------------------------------------------------------------
  !-----     TYPE-BOUND PROCEDURES FOR THE EXPONENTIAL PROPAGATOR     -----
  !------------------------------------------------------------------------

  subroutine direct_solver_lyap(self, vec_in, vec_out)
    !> Linear Operator.
    class(exponential_prop_lyap), intent(in)  :: self
    !> Input vector.
    class(abstract_vector) , intent(in)  :: vec_in
    !> Output vector.
    class(abstract_vector) , intent(out) :: vec_out

    !> Time-integrator.
    type(rks54_class) :: prop
    real(kind=wp)     :: dt = 0.1_wp

    select type(vec_in)
    type is(state_vector_lyap)
       select type(vec_out)
       type is(state_vector_lyap)

          !> Initialize propagator.
          call prop%initialize(n=nx**2, f=rhs_lyap)
          !> Integrate forward in time.
          call prop%integrate(0.0_wp, vec_in%state, dt, self%tau, vec_out%state)

       end select
    end select
    return
  end subroutine direct_solver_lyap

  subroutine adjoint_solver_lyap(self, vec_in, vec_out)
    !> Linear Operator.
    class(exponential_prop_lyap), intent(in)  :: self
    !> Input vector.
    class(abstract_vector) , intent(in)  :: vec_in
    !> Output vector.
    class(abstract_vector) , intent(out) :: vec_out

    !> Time-integrator.
    type(rks54_class) :: prop
    real(kind=wp)     :: dt = 0.1_wp

    select type(vec_in)
    type is(state_vector_lyap)
       select type(vec_out)
       type is(state_vector_lyap)

          !> Initialize propagator.
          call prop%initialize(n=nx**2, f=adjoint_rhs_lyap)
          !> Integrate forward in time.
          call prop%integrate(0.0_wp, vec_in%state, dt, self%tau, vec_out%state)

       end select
    end select
    return
  end subroutine adjoint_solver_lyap

end module Laplacian_Lyapunov
