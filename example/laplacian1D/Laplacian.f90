module Laplacian
  !> RKLIB module for time integration.
  use rklib_module
  !> exmplib module for exponential propagator
  use Lightrom_expmlib
  !> LightKrylov for linear algebra.
  use LightKrylov
  !use LightROM
  !> Standard Library.
  use stdlib_math, only : linspace
  implicit none

  private
  public :: nx, dx, L
  public :: initialize_mesh, laplacian_vec

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

  !-----------------------------------
  !-----     MATVEC LAPLACIAN    -----
  !-----------------------------------

  type, extends(abstract_linop), public :: laplace_operator
   contains
     private
     procedure, pass(self), public :: matvec  => direct_matvec_laplace_operator
     procedure, pass(self), public :: rmatvec => direct_matvec_laplace_operator     ! dummy since Lyapunov equation for Laplacian is symmetric
  end type laplace_operator

  !------------------------------------------
  !-----     EXPONENTIAL PROPAGATOR     -----
  !------------------------------------------

  type, extends(abstract_linop), public :: rklib_exptA_laplacian_vec
     real(kind=wp) :: tau ! Integration time.
   contains
     private
     procedure, pass(self), public :: matvec  => direct_solver_vec
     procedure, pass(self), public :: rmatvec => direct_solver_vec          ! dummy
  end type rklib_exptA_laplacian_vec

  !------------------------------------------------------
  !-----     EXPONENTIAL PROPAGATOR KRYLOV VECTOR   -----
  !------------------------------------------------------

  type, extends(abstract_linop), public :: krylov_exptA_laplacian_vec
     !> operator to be exponentiated
     class(laplace_operator), allocatable :: A
     real(kind=wp) :: tau ! Integration time.
   contains
     private
     procedure, pass(self), public :: matvec  => direct_solver_expm_laplace_vec
     procedure, pass(self), public :: rmatvec => direct_solver_expm_laplace_vec      ! dummy
  end type krylov_exptA_laplacian_vec

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

  !---------------------------
  !-----    Laplacian    -----
  !---------------------------

  subroutine laplacian_vec(vec_out, vec_in)
   
   !> State vector.
   real(kind=wp)  , dimension(:), intent(in)  :: vec_in
   !> Time-derivative.
   real(kind=wp)  , dimension(:), intent(out) :: vec_out

   !> Internal variables.
   integer :: i
   
   !> Left most boundary points 
   vec_out(1) = (- 2*vec_in(1) + vec_in(2)) / dx**2   !> Diffusion term (zero Dirichlet)
   !> Interior nodes.
   do i = 2, nx-1
      vec_out(i) = (vec_in(i+1) - 2*vec_in(i) + vec_in(i-1)) / dx**2
   enddo
   !> Right most boundary points
   vec_out(nx) = (vec_in(nx-1) - 2*vec_in(nx)) / dx**2
      
   return
 end subroutine laplacian_vec

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
    call laplacian_vec(f(1:nx), x(1:nx))
    
    return
  end subroutine rhs

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

  !---------------------------------------------------------------
  !-----     TYPE-BOUND PROCEDURES FOR LAPLACIAN OPERATOR    -----
  !---------------------------------------------------------------

  subroutine direct_matvec_laplace_operator(self, vec_in, vec_out)
   !> Linear Operator.
   class(laplace_operator),intent(in)  :: self
   !> Input vector.
   class(abstract_vector) , intent(in)  :: vec_in
   !> Output vector.
   class(abstract_vector) , intent(out) :: vec_out

   select type(vec_in)
   type is(state_vector)
      select type(vec_out)
      type is(state_vector)
         call laplacian_vec(vec_out%state,vec_in%state)
      end select
   end select
 
   return
 end subroutine direct_matvec_laplace_operator

  !-----------------------------------------------------------------------------
  !-----     TYPE-BOUND PROCEDURES FOR THE EXPONENTIAL PROPAGATOR RKLIB    -----
  !-----------------------------------------------------------------------------

  subroutine direct_solver_vec(self, vec_in, vec_out)
    !> Linear Operator.
    class(rklib_exptA_laplacian_vec), intent(in)  :: self
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
  end subroutine direct_solver_vec

  !--------------------------------------------------------------------------------
 !-----     TYPE-BOUND PROCEDURES FOR THE KRYLOV EXPONENTIAL PROPAGATOR     -----
 !-------------------------------------------------------------------------------

 subroutine direct_solver_expm_laplace_vec(self, vec_in, vec_out)
   !> Linear Operator.
   class(krylov_exptA_laplacian_vec), intent(in)  :: self
   !> Input vector.
   class(abstract_vector) , intent(in)  :: vec_in
   !> Output vector.
   class(abstract_vector) , intent(out) :: vec_out

   !> Internals
   class(state_vector), allocatable  :: tmp1, tmp2
   real(kind=wp), parameter     :: dt0 = 0.1_wp      ! time-step
   real(kind=wp), parameter     :: tol = 1e-9_wp     ! tolerance for the krylov approximation
   real(kind=wp)                :: Tend, dt_last
   integer                      :: i, nsteps, info

   nsteps = floor(self%tau/dt0)
   dt_last = self%tau - nsteps*dt0
   if ( dt_last .lt. atol ) then
      dt_last = dt0
      nsteps = nsteps - 1
   endif
   Tend = 0.0_wp

   select type(vec_in)
   type is(state_vector)
     select type(vec_out)
     type is(state_vector)
       if ( nsteps .eq. 0 ) then          ! nsteps = 0 --> one step (dt_last)
         call kexpm(vec_out, self%A, vec_in, dt_last, tol, info)
         Tend = Tend + dt_last
       elseif ( nsteps .eq. 1 ) then      ! nsteps = 1 --> 2 steps
         allocate(tmp1, source=vec_in)
         call kexpm(tmp1,    self%A, vec_in, dt0,     tol, info)
         Tend = Tend + dt0
         call kexpm(vec_out, self%A, tmp1,   dt_last, tol, info)
         Tend = Tend + dt_last
       else                               ! nsteps > 1
         allocate(tmp1, source=vec_in)
         allocate(tmp2, source=vec_in)
         if ( mod(nsteps,2) .eq.0 ) then  ! nsteps even
           tmp2%state = vec_in%state
           do i = 1, nsteps/2
              call kexpm(tmp1, self%A, tmp2, dt0, tol, info) !> instead of copying data around we 
              call kexpm(tmp2, self%A, tmp1, dt0, tol, info) !  swap input and output
              Tend = Tend + dt0*2.0_wp
           end do
         else                             ! nsteps odd
          tmp1%state = vec_in%state ! inverse tmp1 and tmp2 to end up with output in tmp2
           do i = 1, (nsteps-1)/2
              call kexpm(tmp2, self%A, tmp1, dt0, tol, info) !> instead of copying data around we 
              call kexpm(tmp1, self%A, tmp2, dt0, tol, info) !  swap input and output
              Tend = Tend + dt0*2.0_wp
           end do
           call kexpm(tmp2, self%A, tmp1, dt0, tol, info)
           Tend = Tend + dt0
         endif                            ! odd/even
         !> last step to match Tend and to put output in vec_out
         call kexpm(tmp1, self%A, tmp2, dt_last, tol, info) 
         Tend = Tend + dt_last
         vec_out%state = tmp1%state
       endif                              ! nsteps    
     end select
   end select
     
   return
 end subroutine direct_solver_expm_laplace_vec

end module Laplacian
