module Laplacian2D_LTI_Riccati
   use Laplacian2D_LTI
   !> RKLIB module for time integration.
   use rklib_module
   !> LightKrylov for linear algebra.
   use LightKrylov
   !> Standard Libraries.
   use stdlib_linalg, only: eye
   use stdlib_optval, only: optval
   implicit none

   private
   public :: CARE, laplacian_mat

   !-------------------------------------------
   !-----     LIGHTKRYLOV VECTOR TYPE     -----
   !-------------------------------------------

   type, extends(abstract_vector), public :: state_matrix
      real(kind=wp) :: state(N**2) = 0.0_wp
   contains
      private
      procedure, pass(self), public :: zero
      procedure, pass(self), public :: dot
      procedure, pass(self), public :: scal
      procedure, pass(self), public :: axpby
      procedure, pass(self), public :: rand
   end type state_matrix

   type, extends(abstract_linop), public :: rklib_riccati_mat
      real(kind=wp) :: tau ! Integration time.
   contains
      private
      procedure, pass(self), public :: matvec  => direct_solver_riccati_mat
      procedure, pass(self), public :: rmatvec => direct_solver_riccati_mat  ! dummy, not used
   end type rklib_riccati_mat

contains

   function CARE(X,A,Q,BRinvBT) result(Y)
      real(kind=wp), dimension(n,n) :: X, A, Q, BRinvBT, Y
      Y = matmul(transpose(A), X) + matmul(X, A) + Q - matmul(X, matmul(BRinvBT, X))
   end function CARE

!-----     TYPE-BOUND PROCEDURE FOR VECTORS     -----

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

   subroutine rand(self, ifnorm)
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
   end subroutine rand

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
      real(kind=wp), dimension(N,N) :: mat, dmat
      
      !> Deal with optional argument
      trans = optval(transpose,.false.)
      
      !> Sets the internal variables.
      mat  = reshape(flat_mat_in(1:N**2),(/N, N/))
      dmat = 0.0_wp
      
      if (trans) then
          do j = 1,N
             call laplacian(dmat(j,:), mat(j,:))
          end do
      else
          do j = 1,N
             call laplacian(dmat(:,j), mat(:,j))
          end do
      endif

      !> Reshape for output
      flat_mat_out = reshape(dmat, shape(flat_mat_in))
       
      return
   end subroutine laplacian_mat

   !--------------------------------------------
   !-----    Riccati equation for RKlib    -----
   !--------------------------------------------

   subroutine rhs_riccati(me, t, x, f)
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
      real(kind=wp), dimension(N,N) :: xm, BRinvBTdata
      real(kind=wp), dimension(N**2) :: dv, dvT, CTQcC

      !> Sets the internal variables.
      dv  = 0.0_wp
      dvT = 0.0_wp

      call laplacian_mat(dv,  x, .false.)       ! A @ X
      call laplacian_mat(dvT, x, .true.)        ! ( A @ X.T ).T

      xm = reshape(x, shape(xm))
      f(1:N**2) = dv + dvT + CTQcC - reshape(matmul(xm, matmul(BRinvBTdata,xm)), shape(CTQcC))

      return
   end subroutine rhs_riccati

   !------------------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE EXPONENTIAL PROPAGATOR     -----
   !------------------------------------------------------------------------

   subroutine direct_solver_riccati_mat(self, vec_in, vec_out)
      !> Linear Operator.
      class(rklib_riccati_mat), intent(in)  :: self
      !> Input vector.
      class(abstract_vector) , intent(in)  :: vec_in
      !> Output vector.
      class(abstract_vector) , intent(out) :: vec_out
      !> Time-integrator.
      type(rks54_class) :: prop
      real(kind=wp)     :: dt = 0.0001_wp

      select type(vec_in)
      type is (state_matrix)
         select type(vec_out)
         type is(state_matrix)
            !> Initialize propagator.
            call prop%initialize(n=N**2, f=rhs_riccati)
            !> Integrate forward in time.
            call prop%integrate(0.0_wp, vec_in%state, dt, self%tau, vec_out%state)
         end select
      end select

      return
   end subroutine direct_solver_riccati_mat

end module Laplacian2D_LTI_Riccati