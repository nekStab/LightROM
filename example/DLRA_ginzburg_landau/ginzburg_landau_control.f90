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

   type, extends(rks54_class) :: rks54_class_with_control
      real(dp), allocatable :: K(:,:)
      real(dp), allocatable :: B(:,:)
      real(dp), allocatable :: u_control(:)
      logical :: control_enabled = .false.
      logical, private :: initialised = .false.
   contains
      private
      procedure, pass(self), public :: is_initialised => controller_is_initialised
      procedure, pass(self), public :: setup => setup_controller
      procedure, pass(self), public :: eval => eval_control
   end type rks54_class_with_control

   !---------------------------------------------------------------
   !-----     EXPONENTIAL PROPAGATOR WITH FEEDBACK CONTROL    -----
   !---------------------------------------------------------------

   type, extends(abstract_linop_rdp), public :: exponential_prop_with_control
      real(dp) :: tau ! Integration time.
      type(rks54_class_with_control) :: prop
      logical :: control_enabled = .false.
      logical, private :: initialised = .false.
   contains
      private
      procedure, pass(self), public :: is_initialised => exptA_is_initialised
      procedure, pass(self), public :: init => init_exptA_with_control
      procedure, pass(self), public :: matvec => direct_solver_with_control
      procedure, pass(self), public :: rmatvec => direct_solver_with_control ! dummy
   end type exponential_prop_with_control

contains

   !-------------------------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE exponential_prop_with_control     -----
   !-------------------------------------------------------------------------------

   pure logical function exptA_is_initialised(self) result(initialised)
      class(exponential_prop_with_control), intent(in) :: self
      !! Controller
      initialised = self%initialised
   end function exptA_is_initialised

   subroutine init_exptA_with_control(self, X, B, Rinv, enable_control)
      class(exponential_prop_with_control), intent(inout)  :: self
      !! Linear Operator.
      class(abstract_sym_low_rank_state_rdp), intent(in) :: X
      !! Riccati matrix representation
      class(abstract_vector_rdp), intent(in) :: B(:)
      !! Input matrix
      real(dp), intent(in) :: Rinv(:, :)
      !! Inverse control cost
      logical, optional, intent(in) :: enable_control
      !! enable control?
      ! internals
      character(len=*), parameter :: this_procedure = 'init_exptA_with_control'
      ! setup internal propagator to compute K
      call self%prop%setup(X, B, Rinv)
      ! Choose whether to enable control
      self%control_enabled = optval(enable_control, .true.)
      ! set initialisation flag
      self%initialised = .true.
   end subroutine init_exptA_with_control

   subroutine direct_solver_with_control(self, vec_in, vec_out)
      ! Linear Operator.
      class(exponential_prop_with_control), intent(inout)  :: self
      ! Input vector.
      class(abstract_vector_rdp),  intent(in)  :: vec_in
      ! Output vector.
      class(abstract_vector_rdp),  intent(out) :: vec_out
  
      ! Time-integrator.
      character(len=*), parameter :: this_procedure = 'direct_solver_with_control'
      real(dp) :: dt = 1.0_dp
  
      select type(vec_in)
      type is(state_vector)
         select type(vec_out)
         type is(state_vector)

            ! Initialize propagator.
            call self%prop%initialize(n=2*nx, f=rhs_with_control)
           
            ! Integrate forward in time.
            call self%prop%integrate(0.0_dp, vec_in%state, dt, self%tau, vec_out%state)
      
         class default
               call type_error('vec_out', 'state_vector', 'OUT', this_module, this_procedure)
         end select
      class default
         call type_error('vec_in', 'state_vector', 'IN', this_module, this_procedure)
      end select
   end subroutine direct_solver_with_control

   !--------------------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE rks54_class_with_control     -----
   !--------------------------------------------------------------------------

   pure logical function controller_is_initialised(self) result(initialised)
      class(rks54_class_with_control), intent(in) :: self
      !! Controller
      initialised = self%initialised
   end function controller_is_initialised

   subroutine setup_controller(self, X, B, Rinv)
      class(rks54_class_with_control), intent(inout) :: self
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
            allocate(self%K(size(B), 2*nx), source=zero_rdp)
            allocate(wrk(2*nx, size(B)), source=zero_rdp)
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
      class(rks54_class_with_control), intent(inout) :: self
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
      type is(rks54_class_with_control)
         f = zero_rdp
         call direct_GL(x, f)
         if (me%control_enabled) then
            ! add control
            call me%eval(me%u_control, x)
            f = f + matmul(me%B, me%u_control)
         end if
      class default
         call type_error('me', 'rks54_class_with_control', 'INOUT', this_module, this_procedure)
      end select
   end subroutine rhs_with_control

end module Ginzburg_Landau_Control