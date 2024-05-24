module Laplacian2D_LTI_Riccati_RKlib
   use Laplacian2D_LTI_Riccati_Base
   use laplacian2D_LTI_Riccati_Operators
   !> RKLIB module for time integration.
   use rklib_module
   !> LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov, only : wp => dp
   !> Standard Library.
   use stdlib_math, only : linspace
   use stdlib_optval, only : optval
   use stdlib_linalg, only : eye
   implicit none

   private

   !-----------------------------------------------
   !-----     EXPONENTIAL PROPAGATOR RKLIB    -----
   !-----------------------------------------------

   type, extends(abstract_linop_rdp), public :: rklib_exptA_laplacian
      real(wp) :: tau ! Integration time.
   contains
      private
      procedure, pass(self), public :: matvec  => direct_solver_vec
      procedure, pass(self), public :: rmatvec => direct_solver_vec                  ! dummy
   end type rklib_exptA_laplacian

   type, extends(abstract_linop_rdp), public :: rklib_riccati_mat
      real(wp) :: tau ! Integration time.
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
      class(rk_class), intent(inout)             :: me
      !> Current time.
      real(wp),      intent(in)                :: t
      !> State vector.
      real(wp),        dimension(:), intent(in)  :: x
      !> Time-derivative.
      real(wp),        dimension(:), intent(out) :: f

      f = 0.0_wp
      call laplacian(f(1:N), x(1:N))
      
      return
   end subroutine rhs

   subroutine direct_solver_vec(self, vec_in, vec_out)
      !> Linear Operator.
      class(rklib_exptA_laplacian), intent(in)  :: self
      !> Input vector.
      class(abstract_vector_rdp),   intent(in)  :: vec_in
      !> Output vector.
      class(abstract_vector_rdp),   intent(out) :: vec_out

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

    !-----    Laplacian Riccati RHS for RKlib    -----

   subroutine rhs_riccati(me, t, x, f)
      !> Time-integrator.
      class(rk_class), intent(inout)             :: me
      !> Current time.
      real(wp),        intent(in)                :: t
      !> State vector.
      real(wp),        dimension(:), intent(in)  :: x
      !> Time-derivative.
      real(wp),        dimension(:), intent(out) :: f

      !> Internal variables.
      integer :: i, j, k
      real(wp), dimension(N,N) :: xm
      real(wp), dimension(N**2) :: dv, dvT

      !> Sets the internal variables.
      dv  = 0.0_wp
      dvT = 0.0_wp

      call laplacian_mat(dv,  x, .false.)       ! A @ X
      call laplacian_mat(dvT, x, .true.)        ! ( A @ X.T ).T

      xm = reshape(x, shape(xm))
      f(1:N**2) = dv + dvT + CTQcC - reshape(matmul(xm, matmul(BRinvBTdata,xm)), shape(CTQcC))

      return
   end subroutine rhs_riccati

   subroutine direct_solver_riccati_mat(self, vec_in, vec_out)
      !> Linear Operator.
      class(rklib_riccati_mat),   intent(in)  :: self
      !> Input vector.
      class(abstract_vector_rdp), intent(in)  :: vec_in
      !> Output vector.
      class(abstract_vector_rdp), intent(out) :: vec_out
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

end module Laplacian2D_LTI_Riccati_RKlib