module Laplacian2D_LTI_Riccati
   use Laplacian2D_LTI
   use Laplacian2D_LTI_Lyapunov
   !> RKLIB module for time integration.
   use rklib_module
   !> LightKrylov for linear algebra.
   use LightKrylov
   !> Standard Libraries.
   use stdlib_linalg, only: eye
   use stdlib_optval, only: optval
   implicit none

   private
   public :: rhs_riccati

   type, extends(abstract_linop), public :: rklib_riccati_mat
      real(kind=wp) :: tau ! Integration time.
   contains
      private
      procedure, pass(self), public :: matvec  => direct_solver_riccati_mat
      procedure, pass(self), public :: rmatvec => direct_solver_riccati_mat                ! dummy, not used
   end type rklib_riccati_mat

   !type, extends(abstract_linop), public :: LR_Q
   !contains
   !   private
   !   procedure, pass(self), public :: matvec  => LR_Q_mat
   !   procedure, pass(self), public :: rmatvec => LR_Q_mat                ! dummy since symmetric
   !end type LR_Q

contains

   !function CARE(X,A,Q,BRinvBT) result(Y)
   !   real(kind=wp), dimension(n,n) :: X, A, Q, BRinvBT, Y
   !   Y = matmul(transpose(Adata), X) + matmul(X, Adata) + Qdata - matmul(X, matmul(BRinvBTdata, X))
   !end function CARE

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
      real(kind=wp), dimension(N**2) :: dv, dvT 

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
      real(kind=wp)     :: dt = 0.1_wp

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

   !--------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR LR_Q     -----
   !--------------------------------------------------

   !subroutine LR_Q_mat(self, vec_in, vec_out)
   !   !> Linear Operator.
   !   class(LR_Q) ,            intent(in)  :: self
   !   !> Input vector.
   !   class(abstract_vector) , intent(in)  :: vec_in
   !   !> Output vector.
   !   class(abstract_vector) , intent(out) :: vec_out
   !
   !   ! internal variables
   !   integer :: i
   !   real(kind=wp) , allocatable         :: wrk(:,:)
   !   class(abstract_vector), allocatable :: Uwrk
   !   select type(vec_in)
   !   type is (state_vector)
   !      select type(vec_out)
   !      type is (state_vector)
   !         allocate(wrk(1:size(CT),1))
   !         do i = 1, size(CT)
   !            wrk(i, 1) = CT(i)%dot(vec_in)
   !         end do
   !         allocate (Uwrk, source=CT(1))
   !         call get_vec(Uwrk, CT, matmul(Qcdata, wrk(:,1)))
   !         call vec_out%axpby(0.0_wp, Uwrk, 1.0_wp)
   !      end select
   !   end select
   !
   !   return
   !end subroutine LR_Q_mat


end module Laplacian2D_LTI_Riccati