module Laplacian_Lyapunov
  use Laplacian
  !> RKLIB module for time integration.
  use rklib_module
  !> exmplib module for exponential propagator
  use Lightrom_expmlib
  !> LightKrylov for linear algebra.
  use LightKrylov
  !use LightROM
  !> Standard Libraries.
  use stdlib_math, only : linspace
  use stdlib_optval, only: optval
  implicit none

  private
  public :: initialize_mesh_matrix_eq, laplacian_mat, add_q
  
  !-------------------------------------------
  !-----     LIGHTKRYLOV VECTOR TYPE     -----
  !-------------------------------------------

  type, extends(abstract_vector), public :: state_matrix
     real(kind=wp) :: state(nx**2) = 0.0_wp
   contains
     private
     procedure, pass(self), public :: zero
     procedure, pass(self), public :: dot
     procedure, pass(self), public :: scal
     procedure, pass(self), public :: axpby
  end type state_matrix

  !------------------------------
  !-----     MATVEC LYAP    -----
  !------------------------------

  type, extends(abstract_linop), public :: lyapunov_operator
   contains
     private
     procedure, pass(self), public :: matvec  => direct_matvec_lyapunov_operator
     procedure, pass(self), public :: rmatvec => direct_matvec_lyapunov_operator     ! dummy since Lyapunov equation for Laplacian is symmetric
  end type lyapunov_operator

  !------------------------------------------
  !-----     EXPONENTIAL PROPAGATOR     -----
  !------------------------------------------

  type, extends(abstract_linop), public :: rklib_exptA_lyapunov_mat
     real(kind=wp) :: tau ! Integration time.
   contains
     private
     procedure, pass(self), public :: matvec  => direct_solver_mat
     procedure, pass(self), public :: rmatvec => direct_solver_mat                ! dummy
  end type rklib_exptA_lyapunov_mat

  !------------------------------------------------------
  !-----     EXPONENTIAL PROPAGATOR KRYLOV MATRIX   -----
  !------------------------------------------------------

  type, extends(abstract_linop), public :: krylov_exptA_lyapunov_mat
     !> operator to be exponentiated
     class(lyapunov_operator), allocatable :: A
     real(kind=wp) :: tau ! Integration time.
   contains
     private
     procedure, pass(self), public :: matvec  => direct_solver_expm_lyapunov_mat
     procedure, pass(self), public :: rmatvec => direct_solver_expm_lyapunov_mat      ! dummy
  end type krylov_exptA_lyapunov_mat

contains

  !=========================================
  !=========================================
  !=====                               =====
  !=====     LAPLACIAN LYAPUNOV EQ     =====
  !=====                               =====
  !=========================================
  !=========================================

  !---------------------------------------
  !-----     CONSTRUCT THE MESH      -----
  !---------------------------------------

  subroutine initialize_mesh_matrix_eq()
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
  end subroutine initialize_mesh_matrix_eq

  subroutine add_Q(vec)
   
    !> State vector.
    real(kind=wp), dimension(:), intent(inout)  :: vec
 
    !> Internal variables.
    real(kind=wp), dimension(nx**2) :: q
 
    q = 1.0_wp
    vec = vec + q
    
    return
  end subroutine add_Q

  subroutine laplacian_mat(flat_mat_out, flat_mat_in, transpose)
   
    !> State vector.
    real(kind=wp)  , dimension(:), intent(in)  :: flat_mat_in
    !> Time-derivative.
    real(kind=wp)  , dimension(:), intent(out) :: flat_mat_out
    !> Transpose
    logical, optional :: transpose
    logical           :: trans
 
    !> Internal variables.
    integer :: j
    real(kind=wp), dimension(nx,nx) :: mat, dmat
 
    !> Deal with optional argument
    trans = optval(transpose,.false.)
 
    !> Sets the internal variables.
    mat  = reshape(flat_mat_in(1:nx**2),(/nx, nx/))
    dmat = 0.0_wp
    
    if (trans) then
       do j = 1,nx
          call laplacian_vec(dmat(j,:), mat(j,:))
       end do
    else
       do j = 1,nx
          call laplacian_vec(dmat(:,j), mat(:,j))
       end do
    endif
    
    !> Reshape for output
    flat_mat_out = reshape(dmat, shape(flat_mat_in))
      
   return
  end subroutine laplacian_mat

  !---------------------------------------------------------------
  !-----    Lyapunov equation for the Laplacian for RKlib    -----
  !---------------------------------------------------------------

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
    real(kind=wp), dimension(nx**2) :: dv, dvT, q

    !> Sets the internal variables.
    dv  = 0.0_wp
    dvT = 0.0_wp
    q   = 1.0_wp

    !> We compute the action of the Lyapunov operator without using the transpose of A
    !           L(X) = A @ X +   X @ A.T     + Q
    !                = A @ X + ( A @ X.T ).T + Q
    call laplacian_mat(dv,  x, .false.)       ! A @ X
    call laplacian_mat(dvT, x, .true.)        ! ( A @ X.T ).T

    !> Combine the parts and copy result to the output array.
    f(1:nx**2) = dv + dvT + q
    
    return
  end subroutine rhs_lyap

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
    class(state_matrix), intent(inout) :: self
    self%state = 0.0_wp
    return
  end subroutine zero

  real(kind=wp) function dot(self, vec) result(alpha)
    class(state_matrix)   , intent(in) :: self
    class(abstract_vector), intent(in) :: vec
    select type(vec)
    type is(state_matrix)
       alpha = dot_product(self%state, vec%state)
    end select
    return
  end function dot

  subroutine scal(self, alpha)
    class(state_matrix), intent(inout) :: self
    real(kind=wp)      , intent(in)    :: alpha
    self%state = self%state * alpha
    return
  end subroutine scal

  subroutine axpby(self, alpha, vec, beta)
    class(state_matrix)   , intent(inout) :: self
    class(abstract_vector), intent(in)    :: vec
    real(kind=wp)         , intent(in)    :: alpha, beta
    select type(vec)
    type is(state_matrix)
       self%state = alpha*self%state + beta*vec%state
    end select
    return
  end subroutine axpby

  !--------------------------------------------------------------
  !-----     TYPE-BOUND PROCEDURES FOR LYAPUNOV OPERATOR    -----
  !--------------------------------------------------------------

  subroutine direct_matvec_lyapunov_operator(self, vec_in, vec_out)
   !> Linear Operator.
   class(lyapunov_operator),intent(in)  :: self
   !> Input vector.
   class(abstract_vector) , intent(in)  :: vec_in
   !> Output vector.
   class(abstract_vector) , intent(out) :: vec_out

   !> Internals 
   class(state_matrix), allocatable  :: tmp

   select type(vec_in)
   type is(state_matrix)
      select type(vec_out)
      type is(state_matrix)
         allocate(tmp, source=vec_in) ; tmp%state = 0.0_wp
         call laplacian_mat(tmp%state,     vec_in%state, .false.)
         call laplacian_mat(vec_out%state, vec_in%state, .true.)
         call vec_out%axpby(1.0_wp, tmp, 1.0_wp)
         call add_Q(vec_out%state)              ! this can be more complicated than it is now
      end select
   end select
 
   return
 end subroutine direct_matvec_lyapunov_operator

 !------------------------------------------------------------------------
 !-----     TYPE-BOUND PROCEDURES FOR THE EXPONENTIAL PROPAGATOR     -----
 !------------------------------------------------------------------------

 subroutine direct_solver_mat(self, vec_in, vec_out)
   !> Linear Operator.
   class(rklib_exptA_lyapunov_mat), intent(in)  :: self
   !> Input vector.
   class(abstract_vector) , intent(in)  :: vec_in
   !> Output vector.
   class(abstract_vector) , intent(out) :: vec_out
   !> Time-integrator.
   type(rks54_class) :: prop
   real(kind=wp)     :: dt = 0.1_wp
   
   select type(vec_in)
   type is(state_matrix)
      select type(vec_out)
      type is(state_matrix)
         !> Initialize propagator.
         call prop%initialize(n=nx**2, f=rhs_lyap)
         !> Integrate forward in time.
         call prop%integrate(0.0_wp, vec_in%state, dt, self%tau, vec_out%state)
      end select
   end select
   
   return
 end subroutine direct_solver_mat

 !--------------------------------------------------------------------------------
 !-----     TYPE-BOUND PROCEDURES FOR THE KRYLOV EXPONENTIAL PROPAGATOR     -----
 !-------------------------------------------------------------------------------

 subroutine direct_solver_expm_lyapunov_mat(self, vec_in, vec_out)
   !> Linear Operator.
   class(krylov_exptA_lyapunov_mat), intent(in)  :: self
   !> Input vector.
   class(abstract_vector) , intent(in)  :: vec_in
   !> Output vector.
   class(abstract_vector) , intent(out) :: vec_out

   !> Internals
   class(state_matrix), allocatable  :: tmp1, tmp2
   real(kind=wp), parameter     :: dt0 = 0.01_wp      ! time-step
   real(kind=wp), parameter     :: tol = 1e-11_wp     ! tolerance for the krylov approximation
   logical, parameter           :: verb = .false.
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
   type is(state_matrix)
     select type(vec_out)
     type is(state_matrix)
       if ( nsteps .eq. 0 ) then          ! nsteps = 0 --> one step (dt_last)
         call kexpm(vec_out, self%A, vec_in, dt_last, tol, info, verbosity = verb)
         Tend = Tend + dt_last
         write(*,*) Tend, dt_last, tol, info
       elseif ( nsteps .eq. 1 ) then      ! nsteps = 1 --> 2 steps
         allocate(tmp1, source=vec_in)
         call kexpm(tmp1,    self%A, vec_in, dt0,     tol, info, verbosity = verb)
         Tend = Tend + dt0
         call kexpm(vec_out, self%A, tmp1,   dt_last, tol, info, verbosity = verb)
         Tend = Tend + dt_last
       else                               ! nsteps > 1
         allocate(tmp1, source=vec_in)
         allocate(tmp2, source=vec_in)
         if ( mod(nsteps,2) .eq.0 ) then  ! nsteps even
           tmp2%state = vec_in%state
           do i = 1, nsteps/2
              call kexpm(tmp1, self%A, tmp2, dt0, tol, info, verbosity = verb) !> instead of copying data around we 
              call kexpm(tmp2, self%A, tmp1, dt0, tol, info, verbosity = verb) !  swap input and output
              Tend = Tend + dt0*2.0_wp
           end do
         else                             ! nsteps odd
          tmp1%state = vec_in%state ! inverse tmp1 and tmp2 to end up with output in tmp2
           do i = 1, (nsteps-1)/2
              call kexpm(tmp2, self%A, tmp1, dt0, tol, info, verbosity = verb) !> instead of copying data around we 
              call kexpm(tmp1, self%A, tmp2, dt0, tol, info, verbosity = verb) !  swap input and output
              Tend = Tend + dt0*2.0_wp
           end do
           call kexpm(tmp2, self%A, tmp1, dt0, tol, info, verbosity = verb)
           Tend = Tend + dt0
         endif                            ! odd/even
         !> last step to match Tend and to put output in vec_out
         call kexpm(tmp1, self%A, tmp2, dt_last, tol, info, verbosity = verb) 
         Tend = Tend + dt_last
         vec_out%state = tmp1%state
       endif                              ! nsteps    
     end select
   end select
     
   return
 end subroutine direct_solver_expm_lyapunov_mat

end module Laplacian_Lyapunov
