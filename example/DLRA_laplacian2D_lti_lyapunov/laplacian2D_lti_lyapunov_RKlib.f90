module Laplacian2D_LTI_Lyapunov_RKlib
   !> Standard Library.
   use stdlib_math, only : linspace
   use stdlib_optval, only : optval
   use stdlib_linalg, only : eye
   !> RKLIB module for time integration.
   use rklib_module
   !> LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov, only : wp => dp
   !> Laplacian
   use Laplacian2D_LTI_Lyapunov_Base
   use laplacian2D_LTI_Lyapunov_Operators
   implicit none

   private :: this_module

   character(len=*), parameter :: this_module = 'Laplacian2D_LTI_Lyapunov_RKLib'

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

   type, extends(abstract_linop_rdp), public :: rklib_lyapunov_mat
      real(wp) :: tau ! Integration time.
   contains
      private
      procedure, pass(self), public :: matvec  => direct_solver_mat
      procedure, pass(self), public :: rmatvec => direct_solver_mat                ! dummy
   end type rklib_lyapunov_mat

contains

   !-----------------------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE EXPONENTIAL PROPAGATOR RKLIB    -----
   !-----------------------------------------------------------------------------
   
   !-----    Laplacian RHS for RKlib    -----
   
   subroutine rhs(me, t, x, f)
      !> Time-integrator.
      class(rk_class), intent(inout)      :: me
      !> Current time.
      real(wp), intent(in)                :: t
      !> State vector.
      real(wp), dimension(:), intent(in)  :: x
      !> Time-derivative.
      real(wp), dimension(:), intent(out) :: f

      f = 0.0_wp
      call laplacian(f(1:N), x(1:N))
      
      return
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
      real(wp)          :: dt = 1.0_wp

      select type(vec_in)
      type is (state_vector)
         select type(vec_out)
         type is (state_vector)

             !> Initialize propagator.
             call prop%initialize(n=N, f=rhs)
             !> Integrate forward in time.
             call prop%integrate(0.0_wp, vec_in%state, dt, self%tau, vec_out%state)

         class default
            call stop_error('vec_out must be a state_vector', this_module, 'direct_solver_vec')
         end select
      class default
         call stop_error('vec_in must be a state_vector', this_module, 'direct_solver_vec')
      end select
      return
   end subroutine direct_solver_vec

   !-----    Lyapunov RHS for RKlib    -----

   subroutine rhs_lyap(me, t, x, f)
      !> Time-integrator.
      class(rk_class), intent(inout)      :: me
      !> Current time.
      real(wp), intent(in)                :: t
      !> State vector.
      real(wp), dimension(:), intent(in)  :: x
      !> Time-derivative.
      real(wp), dimension(:), intent(out) :: f

      !> Internal variables.
      integer :: i, j, k
      real(wp), dimension(N**2) :: dv, dvT

      !> Sets the internal variables.
      dv  = 0.0_wp
      dvT = 0.0_wp

      !> We compute the action of the Lyapunov operator without using the transpose of A
      !           L(X) = A @ X +   X @ A.T     + Q
      !                = A @ X + ( A @ X.T ).T + Q
      call laplacian_mat(dv,  x, .false.)       ! A @ X
      call laplacian_mat(dvT, x, .true.)        ! ( A @ X.T ).T

      !> Combine the parts and copy result to the output array.
      f(1:N**2) = dv + dvT + BBT
      
      return
   end subroutine rhs_lyap

   subroutine direct_solver_mat(self, vec_in, vec_out)
      !> Linear Operator.
      class(rklib_lyapunov_mat),  intent(inout)  :: self
      !> Input vector.
      class(abstract_vector_rdp), intent(in)  :: vec_in
      !> Output vector.
      class(abstract_vector_rdp), intent(out) :: vec_out
      !> Time-integrator.
      type(rks54_class) :: prop
      real(wp)          :: dt = 0.1_wp

      select type(vec_in)
      type is (state_matrix)
         select type(vec_out)
         type is(state_matrix)
            !> Initialize propagator.
            call prop%initialize(n=N**2, f=rhs_lyap)
            !> Integrate forward in time.
            call prop%integrate(0.0_wp, vec_in%state, dt, self%tau, vec_out%state)
         class default
            call stop_error('vec_out must be a state_matrix', this_module, 'direct_solver_mat')
         end select
      class default
         call stop_error('vec_in must be a state_matrix', this_module, 'direct_solver_mat')
      end select

      return
   end subroutine direct_solver_mat


end module Laplacian2D_LTI_Lyapunov_RKlib