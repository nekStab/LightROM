module Laplacian2D_LTI_Riccati_RKlib
   !> Standard Library.
   use stdlib_math, only : linspace
   use stdlib_optval, only : optval
   use stdlib_linalg, only : eye
   !> RKLIB module for time integration.
   use rklib_module
   !> LightKrylov for linear algebra.
   use LightKrylov
   !> Laplacian
   use Laplacian2D_LTI_Riccati_Base
   use laplacian2D_LTI_Riccati_Operators   
   implicit none

   private :: this_module
   character(len=*), parameter :: this_module = 'Laplacian2D_LTI_Lyapunov_Base'

   !-----------------------------------------------
   !-----     EXPONENTIAL PROPAGATOR RKLIB    -----
   !-----------------------------------------------

   type, extends(abstract_linop_rdp), public :: rklib_exptA_laplacian
      real(dp) :: tau ! Integration time.
   contains
      private
      procedure, pass(self), public :: matvec  => direct_solver_vec
      procedure, pass(self), public :: rmatvec => direct_solver_vec                  ! dummy
   end type rklib_exptA_laplacian

   type, extends(abstract_linop_rdp), public :: rklib_riccati_mat
      real(dp) :: tau ! Integration time.
   contains
      private
      procedure, pass(self), public :: matvec  => direct_solver_riccati_mat
      procedure, pass(self), public :: rmatvec => direct_solver_riccati_mat  ! dummy, not used
   end type rklib_riccati_mat

contains

   !-----------------------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE EXPONENTIAL PROPAGATOR RKLIB    -----
   !-----------------------------------------------------------------------------
   
   !-----    Laplacian RHS for RKlib    -----
   
   subroutine rhs(me, t, x, f)
      !> Time-integrator.
      class(rk_class), intent(inout)        :: me
      !> Current time.
      real(dp)  , intent(in)                :: t
      !> State vector.
      real(dp)  , dimension(:), intent(in)  :: x
      !> Time-derivative.
      real(dp)  , dimension(:), intent(out) :: f

      f = 0.0_dp
      call laplacian(f(1:N), x(1:N))
      
   end subroutine rhs

   subroutine direct_solver_vec(self, vec_in, vec_out)
      !> Linear Operator.
      class(rklib_exptA_laplacian), intent(inout)  :: self
      !> Input vector.
      class(abstract_vector_rdp),   intent(in)  :: vec_in
      !> Output vector.
      class(abstract_vector_rdp),   intent(out) :: vec_out

      !> Time-integrator.
      type(rks54_class) :: prop
      real(dp)          :: dt = 1.0_dp

      select type(vec_in)
      type is (state_vector)
         select type(vec_out)
         type is (state_vector)

             !> Initialize propagator.
             call prop%initialize(n=N, f=rhs)
             !> Integrate forward in time.
             call prop%integrate(0.0_dp, vec_in%state, dt, self%tau, vec_out%state)

         class default
            call type_error('vec_out', 'state_vector', 'OUT', this_module, this_procedure)
         end select
      class default
         call type_error('vec_in', 'state_vector', 'IN', this_module, this_procedure)
      end select
   end subroutine direct_solver_vec

    !-----    Laplacian Riccati RHS for RKlib    -----

   subroutine rhs_riccati(me, t, x, f)
      !> Time-integrator.
      class(rk_class), intent(inout)        :: me
      !> Current time.
      real(dp)  , intent(in)                :: t
      !> State vector.
      real(dp)  , dimension(:), intent(in)  :: x
      !> Time-derivative.
      real(dp)  , dimension(:), intent(out) :: f

      !> Internal variables.
      integer :: i, j, k
      real(dp), dimension(N,N) :: xm
      real(dp), dimension(N**2) :: dv, dvT

      !> Sets the internal variables.
      dv  = 0.0_dp
      dvT = 0.0_dp

      call laplacian_mat(dv,  x, .false.)       ! A @ X
      call laplacian_mat(dvT, x, .true.)        ! ( A @ X.T ).T

      xm = reshape(x, shape(xm))
      f(1:N**2) = dv + dvT + CTQcC - reshape(matmul(xm, matmul(BRinvBTdata,xm)), [N**2])

   end subroutine rhs_riccati

   subroutine direct_solver_riccati_mat(self, vec_in, vec_out)
      !> Linear Operator.
      class(rklib_riccati_mat),   intent(inout)  :: self
      !> Input vector.
      class(abstract_vector_rdp), intent(in)  :: vec_in
      !> Output vector.
      class(abstract_vector_rdp), intent(out) :: vec_out
      
      !> Time-integrator.
      character(len=*), parameter :: this_procedure = 'direct_solver_mat'
      type(rks54_class) :: prop
      real(dp)          :: dt = 0.1_dp

      select type(vec_in)
      type is (state_matrix)
         select type(vec_out)
         type is(state_matrix)
            !> Initialize propagator.
            call prop%initialize(n=N**2, f=rhs_lyap)
            !> Integrate forward in time.
            call prop%integrate(0.0_dp, vec_in%state, dt, self%tau, vec_out%state)
         class default
            call type_error('vec_out', 'state_matrix', 'OUT', this_module, this_procedure)
         end select
      class default
         call type_error('vec_in', 'state_matrix', 'IN', this_module, this_procedure)
      end select
   end subroutine direct_solver_riccati_mat

end module Laplacian2D_LTI_Riccati_RKlib