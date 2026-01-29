module Ginzburg_Landau_Control
   ! Standard Library.
   use stdlib_optval, only : optval
   ! RKLIB module for time integration.
   use rklib_module
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov_utils, only : assert_shape
   use LightKrylov_Constants, only : zero_rdp
   ! LightROM
   use LightROM_AbstractLTIsystems ! abstract_lti_system
   use LightROM_Utils
   ! Ginzburg Landau
   use Ginzburg_Landau_Base
   use Ginzburg_Landau_Operators
   use Ginzburg_Landau_Utils, only : get_state
   implicit none

   private :: this_module
   character(len=*), parameter :: this_module = 'Ginzburg_Landau_Control'

   !----------------------------------------------------------
   !-----     RK_LIB rks54_class with control inputs     -----
   !----------------------------------------------------------

    type, extends(rks54_class) :: controlled_rks54_class
        real(dp), allocatable :: K(:,:)
        real(dp), allocatable :: B(:,:)
        real(dp), allocatable :: u_control(:)
        logical :: control_enabled = .false.
        logical, private :: initialised = .false.
    contains
        private
        procedure, pass(self), public :: is_initialised
        procedure, pass(self), public :: setup => setup_controller
        procedure, pass(self), public :: eval => eval_control
    end type controlled_rks54_class

    !---------------------------------------------------------------
    !-----     EXPONENTIAL PROPAGATOR WITH FEEDBACK CONTROL    -----
    !---------------------------------------------------------------

    type, extends(abstract_linop_rdp), public :: feedback_exponential_prop
        real(dp) :: tau ! Integration time.
    contains
        private
        procedure, pass(self), public :: matvec => direct_solver_with_control
        procedure, pass(self), public :: rmatvec => direct_solver_with_control ! dummy
    end type feedback_exponential_prop

contains

    pure logical function is_initialised(self) result(initialised)
        class(controlled_rks54_class), intent(in) :: self
        !! Controller
        initialised = self%initialised
    end function

    subroutine setup_controller(self, X, B, Rinv)
        class(controlled_rks54_class), intent(inout) :: self
        !! Controller
        class(abstract_sym_low_rank_state_rdp), intent(in) :: X
        !! Riccati matrix representation
        class(abstract_vector_rdp), intent(in) :: B(:)
        !! Input matrix
        real(dp), intent(in) :: Rinv(:, :)
        ! internals
        character(len=*), parameter :: this_procedure = 'setup_controller'
        class(abstract_vector_rdp), allocatable :: KT(:)
        real(dp), allocatable :: wrk(:,:)
        select type (X)
        type is (LR_state)
            select type (B)
            type is (state_vector)
                ! copy input matrix
                allocate(self%B(2*nx, size(B)), source=zero_rdp)
                call get_state(self%B, B, this_procedure)
                ! compute LQR gains and store them
                call LQR_gain(KT, X, B, Rinv) ! will allocate KT as a state vector
                allocate(self%K(X%rk, 2*nx), source=zero_rdp)
                allocate(wrk(2*nx, X%rk), source=zero_rdp)
                call get_state(wrk, KT, this_procedure)
                self%K = transpose(wrk)
                ! set init flag
                allocate(self%u_control(size(B)))
                self%initialised = .true.
            class default
                call type_error('B', 'state_vector', this_module, this_procedure)
            end select
        class default
            call type_error('X', 'LR_state', this_module, this_procedure)
        end select
    end subroutine setup_controller

    subroutine eval_control(self, u, x)
        !! Evaluate the LQR control law on a state
        class(controlled_rks54_class), intent(inout) :: self
        !! Controller
        real(dp), intent(out) :: u(:)
        !! Control action
        real(dp), intent(in)  :: x(:)
        !! Current state
        ! internals
        character(len=*), parameter :: this_procedure = 'eval_controller'
        if (.not. self%is_initialised()) then
           call stop_error("The controller is not initialised", this_module, this_procedure)
        end if
        u = -matmul(self%K, x)
    end subroutine eval_control

    subroutine rhs_with_control(me, t, x, f)
        class(rk_class), intent(inout) :: me
        real(dp),        intent(in)    :: t
        real(dp),        intent(in)    :: x(:)
        real(dp),        intent(out)   :: f(:)
        ! internala
        character(len=*), parameter :: this_procedure = 'rhs_with_control'
        real(dp), allocatable :: f_control(:)
        select type(me)
        type is(controlled_rks54_class)
            f = zero_rdp
            if (me%control_enabled) then
                call direct_GL(x, f)
                ! add control
                call me%eval(me%u_control, x)
                f = f + matmul(me%B, me%u_control)
            else
                call direct_GL(x, f)
            end if
        class default
            call type_error('me', 'controlled_rks54_class', 'INOUT', this_module, this_procedure)
        end select
    end subroutine rhs_with_control

    subroutine direct_solver_with_control(self, vec_in, vec_out)
        ! Linear Operator.
        class(feedback_exponential_prop),     intent(inout)  :: self
        ! Input vector.
        class(abstract_vector_rdp),  intent(in)  :: vec_in
        ! Output vector.
        class(abstract_vector_rdp),  intent(out) :: vec_out
  
        ! Time-integrator.
        character(len=*), parameter :: this_procedure = 'direct_solver_with_control'
        type(controlled_rks54_class) :: prop
        real(dp)               :: dt = 1.0_dp
  
        select type(vec_in)
        type is(state_vector)
            select type(vec_out)
            type is(state_vector)
        
                ! Initialize propagator.
                call prop%initialize(n=2*nx, f=rhs_with_control)
             
                ! Integrate forward in time.
                call prop%integrate(0.0_dp, vec_in%state, dt, self%tau, vec_out%state)
        
            class default
                call type_error('vec_out', 'state_vector', 'OUT', this_module, this_procedure)
            end select
        class default
            call type_error('vec_in', 'state_vector', 'IN', this_module, this_procedure)
        end select
     end subroutine direct_solver_with_control
     

end module Ginzburg_Landau_Control