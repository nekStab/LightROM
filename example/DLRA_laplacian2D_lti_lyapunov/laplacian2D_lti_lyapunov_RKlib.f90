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

   character*128, parameter :: this_module = 'Laplacian2D_LTI_Lyapunov_RKLib'
   ! exptA
   public :: exptA_rklib
   
   !-----------------------------------------------
   !-----     LIGHTKRYLOV LTI SYSTEM TYPE     -----
   !-----------------------------------------------

   type, extends(abstract_lti_system_rdp), public :: lti_system
   contains
      private
      procedure, pass(self), public :: initialize => initialize_lti_system
   end type lti_system

   !-----------------------------------------------
   !-----     EXPONENTIAL PROPAGATOR RKLIB    -----
   !-----------------------------------------------

   type, extends(abstract_linop_rdp), public :: rklib_exptA_laplacian
      !! Exponential propagator for the Laplacian
      real(wp) :: tau ! Integration time.
   contains
      private
      procedure, pass(self), public :: matvec  => direct_solver_vec
      procedure, pass(self), public :: rmatvec => direct_solver_vec                  ! dummy
   end type rklib_exptA_laplacian

   type, extends(abstract_linop_rdp), public :: rklib_lyapunov_mat
      !! Direct integrator for the Lyapunov equation for the Laplacian
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
      class(rklib_exptA_laplacian), intent(in)  :: self
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

         end select
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
      class(rklib_lyapunov_mat),  intent(in)  :: self
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
         end select
      end select

      return
   end subroutine direct_solver_mat

   !--------------------------------------
   !-----     EXP(tA) SUBROUTINE     -----
   !--------------------------------------

   subroutine exptA_rklib(vec_out, A, vec_in, tau, info, trans)
      !! Subroutine for the exponential propagator that conforms with the abstract interface
      !! defined in expmlib.f90
      class(abstract_vector_rdp),  intent(out)   :: vec_out
      !! Output vector
      class(abstract_linop_rdp),   intent(inout) :: A
      !! Linear operator
      class(abstract_vector_rdp),  intent(in)    :: vec_in
      !! Input vector.
      real(wp),                    intent(in)    :: tau
      !! Integration horizon
      integer,                     intent(out)   :: info
      !! Information flag
      logical, optional,           intent(in)    :: trans
      logical                                    :: transpose
      !! Direct or Adjoint?

      ! optional argument
      transpose = optval(trans, .false.)

      ! time integrator
      select type (vec_in)
      type is (state_vector)
         select type (vec_out)
         type is (state_vector)
            select type (A)
            type is (rklib_exptA_laplacian)
               ! set integration time
               A%tau = tau
               if (transpose) then
                  call A%rmatvec(vec_in, vec_out) ! this is not really necessary since operator is self-adjoint
               else
                  call A%matvec(vec_in, vec_out)
               end if
            end select
         end select
      end select

   end subroutine exptA_rklib


   !--------------------------------------------------------
   !-----     TYPE BOUND PROCEDURES FOR LTI SYSTEMS    -----
   !--------------------------------------------------------

   subroutine initialize_lti_system(self, A, prop, B, CT, D)
      class(lti_system),           intent(inout) :: self
      class(abstract_linop_rdp),   intent(in)    :: A
      class(abstract_linop_rdp),   intent(in)    :: prop
      class(abstract_vector_rdp),  intent(in)    :: B(:)
      class(abstract_vector_rdp),  intent(in)    :: CT(:)
      real(wp),          optional, intent(in)    :: D(:,:)

      ! internal variables
      integer                                :: rk_b, rk_c

      ! Operator
      select type (A)
      type is (laplace_operator)
         allocate(self%A, source=A)
      end select
      ! Exp Prop
      select type (prop)
      type is (rklib_exptA_laplacian)
         allocate(self%prop, source=prop)
      end select
      ! Input
      select type (B)
      type is (state_vector)
         rk_b = size(B)
         allocate(self%B(1:rk_b), source=B(1:rk_b))
      end select
      ! Output
      select type (CT)
         type is (state_vector)
         rk_c = size(CT)
         allocate(self%CT(1:rk_c), source=CT(1:rk_c))
      end select
      ! Throughput
      allocate(self%D(1:rk_c, 1:rk_b))
      if (present(D)) then
         call assert_shape(D, (/ rk_c, rk_b /), 'initialize_lti_system', 'D')
         self%D = D
      else
         self%D = 0.0_wp
      end if
      return
   end subroutine initialize_lti_system


end module Laplacian2D_LTI_Lyapunov_RKlib
