module Laplacian_Lyapunov_kexpm
  use Laplacian
  use Laplacian_Lyapunov
  !> exmplib module for exponential propagator
  use Lightrom_expmlib
  !> LightKrylov for linear algebra.
  use LightKrylov
  !use LightROM
  !> Standard Library.
  use stdlib_math, only : linspace
  implicit none

  private
  public :: initialize_mesh_lyap


  !-------------------------
  !-----     MATVEC    -----
  !-------------------------

  type, extends(abstract_linop), public :: matvec_laplacian
   contains
     private
     procedure, pass(self), public :: matvec  => direct_matvec
     procedure, pass(self), public :: rmatvec => direct_matvec     ! dummy since Laplacian is symmetric
  end type matvec_laplacian

  !------------------------------
  !-----     MATVEC LYAP    -----
  !------------------------------

  type, extends(abstract_linop), public :: matvec_lyap
     class(matvec_laplacian), allocatable :: A
   contains
     private
     procedure, pass(self), public :: matvec  => direct_matvec_lyap
     procedure, pass(self), public :: rmatvec => direct_matvec_lyap     ! dummy since Lyapunov equation for Laplacian is symmetric
  end type matvec_lyap

  !------------------------------------------
  !-----     EXPONENTIAL PROPAGATOR     -----
  !------------------------------------------

  type, extends(abstract_linop), public :: krylov_exponential_prop_lyap
     !> operator to be exponentiated
     class(matvec_lyap), allocatable :: A
     real(kind=wp) :: tau ! Integration time.
   contains
     private
     procedure, pass(self), public :: matvec  => direct_solver_expm_lyap
     procedure, pass(self), public :: rmatvec => direct_solver_expm_lyap      ! dummy
  end type krylov_exponential_prop_lyap

contains

  !=====================================
  !=====================================
  !=====                           =====
  !=====     LAPLACIAN LYAPUNOV    =====
  !=====                           =====
  !=====================================
  !=====================================

  !------------------------------
  !-----      Laplacian     -----
  !------------------------------

  subroutine apply_A(u_out, u_in, trans)
   
    !> State vector.
    real(kind=wp)  , dimension(:), intent(in)  :: u_in
    !> Time-derivative.
    real(kind=wp)  , dimension(:), intent(out) :: u_out
    !> Transposed
    logical ,                      intent(in)  :: trans

    !> Internal variables.
    integer :: i, j
    real(kind=wp), dimension(nx,nx) :: v, dv

    !> Sets the internal variables.
    v   = reshape(u_in(1:nx**2),(/nx, nx/))
    dv  = 0.0_wp

    if (trans) then
      do j = 1,nx
        !> Left most boundary points 
        dv(j,1) = (- 2*v(j,1) + v(j,2)) / dx**2   !> Diffusion term (zero Dirichlet)
        !> Interior nodes.
        do i = 2, nx-1
           dv(j,i) = (v(j,i+1) - 2*v(j,i) + v(j,i-1)) / dx**2
        enddo
        !> Right most boundary points
        dv(j,nx) = (v(j,nx-1) - 2*v(j,nx)) / dx**2
      enddo
    else
      do j = 1,nx
        !> Left most boundary points 
        dv( 1,j) = (- 2*v(1,j) + v(2,j)) / dx**2   !> Diffusion term (zero Dirichlet)
        !> Interior nodes.
        do i = 2, nx-1
           dv( i,j) = (v(i+1,j) - 2*v(i,j) + v(i-1,j)) / dx**2
        enddo
        !> Right most boundary points
        dv( nx,j) = (v(nx-1,j) - 2*v(nx,j)) / dx**2
      enddo
    endif
    
    !> Combine the parts and copy result to the output array.
    u_out = reshape(dv, shape(u_in))
    
    return
  end subroutine apply_A

  subroutine add_Q(vec)
   
   !> State vector.
   real(kind=wp)  , dimension(:), intent(inout)  :: vec

   !> Internal variables.
   real(kind=wp), dimension(nx*2) :: q

   q = 1.0_wp
   vec = vec + q
   
   return
  end subroutine add_Q

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

  !-------------------------------------------------------------
  !-----     TYPE-BOUND PROCEDURES FOR MATVEC LAPLACIAN    -----
  !-------------------------------------------------------------

  subroutine direct_matvec(self, vec_in, vec_out)
   !> Linear Operator.
   class(matvec_laplacian),   intent(in)  :: self
   !> Input vector.
   class(abstract_vector) , intent(in)  :: vec_in
   !> Output vector.
   class(abstract_vector) , intent(out) :: vec_out

   select type(vec_in)
   type is(state_vector_lyap)
      select type(vec_out)
      type is(state_vector_lyap)
         call apply_A(vec_out%state, vec_in%state, .false.)
      end select
   end select
  
   return
 end subroutine direct_matvec

  !--------------------------------------------------------
  !-----     TYPE-BOUND PROCEDURES FOR MATVEC LYAP    -----
  !--------------------------------------------------------

  subroutine direct_matvec_lyap(self, vec_in, vec_out)
   !> Linear Operator.
   class(matvec_lyap),     intent(in)  :: self
   !> Input vector.
   class(abstract_vector) , intent(in)  :: vec_in
   !> Output vector.
   class(abstract_vector) , intent(out) :: vec_out

   !> Internals 
   class(state_vector_lyap), allocatable  :: tmp

   select type(vec_in)
   type is(state_vector_lyap)
      select type(vec_out)
      type is(state_vector_lyap)
         allocate(tmp, source=vec_in) ; tmp%state = 0.0_wp
         call apply_A(tmp%state,     vec_in%state, .false.)
         call apply_A(vec_out%state, vec_in%state, .true. )
         call vec_out%axpby(1.0_wp, tmp, 1.0_wp)
         call add_Q(vec_out%state)              ! this can be more complicated than it is now
      end select
   end select
 
   return
 end subroutine direct_matvec_lyap

  !------------------------------------------------------------------------
  !-----     TYPE-BOUND PROCEDURES FOR THE EXPONENTIAL PROPAGATOR     -----
  !------------------------------------------------------------------------

  subroutine direct_solver_expm_lyap(self, vec_in, vec_out)
    !> Linear Operator.
    class(krylov_exponential_prop_lyap), intent(in)  :: self
    !> Input vector.
    class(abstract_vector) , intent(in)  :: vec_in
    !> Output vector.
    class(abstract_vector) , intent(out) :: vec_out

    !> Internals
    class(state_vector_lyap), allocatable  :: tmp1(:), tmp2(:)
    real(kind=wp), parameter     :: dt0 = 0.1_wp      ! time-step
    real(kind=wp), parameter     :: tol = 1e-9_wp     ! tolerance for the krylov approximation
    real(kind=wp)                :: Tend, dt_last
    integer                      :: i, nsteps, info

    nsteps = floor(self%tau/dt0)
    dt_last = self%tau - nsteps*dt0

    Tend = 0.0_wp
    select type(vec_in)
    type is(state_vector_lyap)
      select type(vec_out)
      type is(state_vector_lyap)
        allocate(tmp1(1), source=vec_in)
        allocate(tmp2(1), source=vec_in)
        if ( mod(nsteps,2) .eq.0 ) then
          tmp2(1)%state = vec_in%state
          do i = 1, nsteps/2
             call kexpm(tmp1, self%A, tmp2, dt0, tol, info) !> instead of copying data around we 
             call kexpm(tmp2, self%A, tmp1, dt0, tol, info) !  swap input and output
             Tend = Tend + dt0*2.0_wp
          end do
        else ! nsteps uneven
         tmp1(1)%state = vec_in%state ! inverse tmp1 and tmp2 to end up with output in tmp2
          do i = 1, (nsteps-1)/2
             call kexpm(tmp2, self%A, tmp1, dt0, tol, info) !> instead of copying data around we 
             call kexpm(tmp1, self%A, tmp2, dt0, tol, info) !  swap input and output
             Tend = Tend + dt0*2.0_wp
          end do
          call kexpm(tmp2, self%A, tmp1, dt0, tol, info)
          Tend = Tend + dt0
        endif
        !> last step to match Tend and to put output in vec_out
        call kexpm(tmp1, self%A, tmp2, dt_last, tol, info) 
        Tend = Tend + dt_last
        vec_out%state = tmp1(1)%state
      end select
    end select
      
    return
  end subroutine direct_solver_expm_lyap

end module Laplacian_Lyapunov_kexpm
