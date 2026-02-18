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

   !-----------------------------------------------------------------
   !-----     RK_LIB rks54_class with control/sensor inputs     -----
   !-----------------------------------------------------------------

   type, extends(rks54_class) :: rks54_class_with_control
      real(dp), allocatable :: gain(:,:)
      real(dp), allocatable :: input(:,:)
      real(dp), allocatable :: u_control(:)
      logical :: adjoint = .false.
      logical :: control_enabled = .false.
      logical, private :: initialised = .false.
   contains
      private
      procedure, pass(self), public :: is_initialised => controller_is_initialised
      procedure, pass(self), public :: check_control_enabled
      procedure, pass(self), public :: setup => setup_controller
      procedure, pass(self), public :: eval => eval_control
   end type rks54_class_with_control

   !---------------------------------------------------------------
   !-----     EXPONENTIAL PROPAGATOR WITH FEEDBACK CONTROL    -----
   !---------------------------------------------------------------

   type, extends(abstract_expta_linop_rdp), public :: exponential_prop_with_control
      !real(dp) :: tau ! Integration time.
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

   !---------------------------------------------------------------------------
   !-----     RK_LIB rks54_class with feedback control with estimation    -----
   !---------------------------------------------------------------------------

   type, extends(rks54_class) :: rks54_class_LQG
      real(dp), allocatable :: K(:,:)
      real(dp), allocatable :: L(:,:)
      real(dp), allocatable :: B(:,:)
      real(dp), allocatable :: C(:,:)
      real(dp), allocatable :: xu_control(:)
      real(dp), allocatable :: xe_control(:)
      real(dp), allocatable :: e_measure(:)
      logical :: control_enabled = .false.
      logical, private :: initialised = .false.
   contains
      private
      procedure, pass(self), public :: is_initialised => LQG_is_initialised
      procedure, pass(self), public :: check_LQG_enabled
      procedure, pass(self), public :: setup => setup_LQG
      procedure, pass(self), public :: eval => eval_LQG_control
   end type rks54_class_LQG

   !------------------------------------------------------------------------------
   !-----     EXPONENTIAL PROPAGATOR WITH FEEDBACK CONTROL AND ESTIMATION    -----
   !------------------------------------------------------------------------------

   type, extends(abstract_exptA_linop_rdp), public :: exponential_prop_LQG
      !real(dp) :: tau ! Integration time.
      type(rks54_class_LQG) :: prop
      logical :: control_enabled = .false.
      logical, private :: initialised = .false.
   contains
      private
      procedure, pass(self), public :: is_initialised => exptA_LQG_is_initialised
      procedure, pass(self), public :: init => init_exptA_LQG
      procedure, pass(self), public :: matvec => direct_solver_LQG
      procedure, pass(self), public :: rmatvec => direct_solver_LQG ! dummy
   end type exponential_prop_LQG

contains

   !--------------------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE rks54_class_with_control     -----
   !--------------------------------------------------------------------------

   pure logical function controller_is_initialised(self) result(initialised)
      class(rks54_class_with_control), intent(in) :: self
      !! Controller
      initialised = self%initialised
   end function controller_is_initialised
   
   subroutine check_control_enabled(self, linop_setting)
      class(rks54_class_with_control), intent(inout) :: self
      !! Controller
      logical, intent(in) :: linop_setting
      !! Setting at the linop level for control enabled/disabled
      self%control_enabled = linop_setting
   end subroutine check_control_enabled

   subroutine setup_controller(self, X, input, W, adjoint)
      class(rks54_class_with_control), intent(inout) :: self
      !! Controller
      class(abstract_sym_low_rank_state_rdp), intent(in) :: X
      !! Riccati matrix representation
      class(abstract_vector_rdp), intent(in) :: input(:)
      !! Input matrix (B or CT)
      real(dp), intent(in) :: W(:, :)
      !! Weights (Inverse control cost Rinv or inverse sensor noise Vinv)
      logical, intent(in) :: adjoint
      !! adjoint?
      ! internals
      character(len=*), parameter :: this_procedure = 'setup_controller'
      class(abstract_vector_rdp), allocatable :: gain(:)
      real(dp), allocatable :: wrk(:,:)
      select type (X)
      type is (LR_state)
         select type (input)
         type is (state_vector)
            ! copy input matrix, compute LQR/LQE gains and store them
            allocate(self%input(N, size(input)), source=zero_rdp)
            allocate(wrk(N, size(input)), source=zero_rdp)
            if (adjoint) then
               ! self%input = C, input = C.T
               call get_state(wrk, input, this_procedure)
               self%input = transpose(wrk)
               ! self%gain = L
               allocate(self%gain(N, size(input)), source=zero_rdp)
               call LQE_gain(gain, X, input, W) ! will allocate gain as a state vector
               call get_state(self%gain, gain, this_procedure)
            else
               ! self%input = B
               call get_state(self%input, input, this_procedure)
               ! self%gain = K, gain = K.T
               allocate(self%gain(size(input), N), source=zero_rdp)
               call LQR_gain(gain, X, input, W) ! will allocate gain as a state vector
               call get_state(wrk, gain, this_procedure)
               self%gain = transpose(wrk)
            end if
            ! set init flag
            allocate(self%u_control(size(input)))
            self%initialised = .true.
            self%adjoint = adjoint
         class default
            call type_error('input', 'state_vector', this_module, this_procedure)
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
      if (self%adjoint) then
         ! u = -Cx
         u = -matmul(self%input, x)
      else
         ! u = -Kx
         u = -matmul(self%gain, x)
      end if
   end subroutine eval_control

   subroutine rhs_with_control(me, t, x, f)
      class(rk_class), intent(inout) :: me
      real(dp),        intent(in)    :: t
      real(dp),        intent(in)    :: x(:)
      real(dp),        intent(out)   :: f(:)
      ! internals
      character(len=*), parameter :: this_procedure = 'rhs_with_control'
      select type(me)
      type is(rks54_class_with_control)
         f = zero_rdp
         call direct_GL(x, f)
         if (me%control_enabled) then
            ! add control
            call me%eval(me%u_control, x)
            if (me%adjoint) then
               ! f = -Lu
               f = f + matmul(me%gain, me%u_control)
            else
               ! f = -Bu
               f = f + matmul(me%input, me%u_control)
            end if
         end if
      class default
         call type_error('me', 'rks54_class_with_control', 'INOUT', this_module, this_procedure)
      end select
   end subroutine rhs_with_control

   !-------------------------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE exponential_prop_with_control     -----
   !-------------------------------------------------------------------------------

   pure logical function exptA_is_initialised(self) result(initialised)
      class(exponential_prop_with_control), intent(in) :: self
      !! Controller
      initialised = self%initialised
   end function exptA_is_initialised

   subroutine init_exptA_with_control(self, X, input, W, adjoint, enable_control)
      class(exponential_prop_with_control), intent(inout)  :: self
      !! Linear Operator.
      class(abstract_sym_low_rank_state_rdp), intent(in) :: X
      !! Riccati matrix representation
      class(abstract_vector_rdp), intent(in) :: input(:)
      !! Input matrix (either B or CT)
      real(dp), intent(in) :: W(:, :)
      !! Weights (Inverse control cost Rinv or inverse sensor noise Vinv)
      logical, optional, intent(in) :: adjoint
      !! adjoint?
      logical, optional, intent(in) :: enable_control
      !! enable control?
      ! internals
      character(len=*), parameter :: this_procedure = 'init_exptA_with_control'
      logical :: adjoint_
      adjoint_ = optval(adjoint, .false.)
      ! Initialize propagator.
      call self%prop%initialize(n=2*nx, f=rhs_with_control)
      ! setup internal propagator to compute K/L
      call self%prop%setup(X, input, W, adjoint)
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
      
            ! check if control should be enabled
            call self%prop%check_control_enabled(self%control_enabled)
            
            ! Integrate forward in time.
            call self%prop%integrate(0.0_dp, vec_in%state, dt, self%tau, vec_out%state)
      
         class default
               call type_error('vec_out', 'state_vector', 'OUT', this_module, this_procedure)
         end select
      class default
         call type_error('vec_in', 'state_vector', 'IN', this_module, this_procedure)
      end select
   end subroutine direct_solver_with_control

   !-----------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE rks54_class_LQG     -----
   !-----------------------------------------------------------------

   pure logical function LQG_is_initialised(self) result(initialised)
      class(rks54_class_LQG), intent(in) :: self
      !! Controller
      initialised = self%initialised
   end function LQG_is_initialised
   
   subroutine check_LQG_enabled(self, linop_setting)
      class(rks54_class_LQG), intent(inout) :: self
      !! Controller
      logical, intent(in) :: linop_setting
      !! Setting at the linop level for control enabled/disabled
      self%control_enabled = linop_setting
   end subroutine check_LQG_enabled

   subroutine setup_LQG(self, Xdir, Xadj, B, CT, Rinv, Vinv)
      class(rks54_class_LQG), intent(inout) :: self
      !! Controller
      class(abstract_sym_low_rank_state_rdp), intent(in) :: Xdir
      !! Riccati matrix representation of solution to direct problem (LQR)
      class(abstract_sym_low_rank_state_rdp), intent(in) :: Xadj
      !! Riccati matrix representation of solution to adjoint problem (LQE)
      class(abstract_vector_rdp), intent(in) :: B(:)
      !! Input matrix
      class(abstract_vector_rdp), intent(in) :: CT(:)
      !! Sensor matrix
      real(dp), intent(in) :: Rinv(:, :)
      !! Inverse control cost
      real(dp), intent(in) :: Vinv(:, :)
      !! Inverse sensor noise
      ! internals
      character(len=*), parameter :: this_procedure = 'setup_LQG'
      class(abstract_vector_rdp), allocatable :: gain(:)
      real(dp), allocatable :: wrk(:,:)
      select type (Xdir)
      type is (LR_state)
         select type (Xadj)
         type is (LR_state)
            select type (B)
            type is (state_vector)
               select type (CT)
               type is (state_vector)
                  ! copy input matrices
                  allocate(self%B(N, size(B)), source=zero_rdp)
                  call get_state(self%B, B, this_procedure)
                  allocate(self%C(N, size(CT)), source=zero_rdp)
                  allocate(wrk(N, size(CT)), source=zero_rdp)
                  call get_state(wrk, CT, this_procedure)
                  self%C = transpose(wrk)
                  deallocate(wrk)
                  ! compute LQE gains and store them
                  call LQE_gain(gain, Xadj, CT, Vinv) ! will allocate gain as a state vector
                  allocate(self%L(N, size(CT)), source=zero_rdp)
                  call get_state(self%L, gain, this_procedure)
                  ! compute LQR gains and store them
                  call LQR_gain(gain, Xdir, B, Rinv) ! will allocate gain as a state vector
                  allocate(self%K(size(B), N), source=zero_rdp)
                  allocate(wrk(N, size(B)), source=zero_rdp)
                  call get_state(wrk, gain, this_procedure)
                  self%K = transpose(wrk)
                  deallocate(wrk)
                  ! set init flag
                  allocate(self%xu_control(size(B)))
                  allocate(self%xe_control(size(B)))
                  allocate(self%e_measure(size(CT)))
                  self%initialised = .true.
               class default
                  call type_error('CT', 'state_vector', this_module, this_procedure)
               end select
            class default
               call type_error('B', 'state_vector', this_module, this_procedure)
            end select
         class default
            call type_error('Xadj', 'LR_state', this_module, this_procedure)
         end select
      class default
         call type_error('Xdir', 'LR_state', this_module, this_procedure)
      end select
   end subroutine setup_LQG

   subroutine eval_LQG_control(self, xu, xe, y, x, e)
      !! Evaluate the LQR control law on a state
      class(rks54_class_LQG), intent(inout) :: self
      !! Controller
      real(dp), intent(out) :: xu(:)
      !! Control action: -Kx
      real(dp), intent(out) :: xe(:)
      !! Feedback on error: Ke
      real(dp), intent(out) :: y(:)
      !! Error measurement: Ce
      real(dp), intent(in)  :: x(:)
      !! Current state
      real(dp), intent(in)  :: e(:)
      !! Current error estimate
      ! internals
      character(len=*), parameter :: this_procedure = 'eval_LQG'
      if (.not. self%is_initialised()) then
         call stop_error("The controller is not initialised", this_module, this_procedure)
      end if
      xu = -matmul(self%K, x)
      xe =  matmul(self%K, e)
      y  = -matmul(self%C, e)
   end subroutine eval_LQG_control

   subroutine rhs_LQG(me, t, x, f)
      class(rk_class), intent(inout) :: me
      real(dp),        intent(in)    :: t
      real(dp),        intent(in)    :: x(:)
      real(dp),        intent(out)   :: f(:)
      ! internals
      character(len=*), parameter :: this_procedure = 'rhs_LQG'
      real(dp) :: state(N), error(N)
      real(dp) :: xdot(N),  edot(N)
      select type(me)
      type is(rks54_class_LQG)
         f = zero_rdp
         ! split imputs
         state = x(:N)
         error = x(N+1:)
         ! action of A
         call direct_GL(state, xdot)
         call direct_GL(error, edot)
         if (me%control_enabled) then
            ! add feedback
            call me%eval(me%xu_control, me%xe_control, me%e_measure, state, error)
            xdot = xdot + matmul(me%B, me%xu_control + me%xe_control)
            edot = edot + matmul(me%L, me%e_measure)
         end if
         f(:N)   = xdot
         f(N+1:) = edot
      class default
         call type_error('me', 'rks54_class_LQG', 'INOUT', this_module, this_procedure)
      end select
   end subroutine rhs_LQG

   !----------------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE exponential_prop_LQG     -----
   !----------------------------------------------------------------------

   pure logical function exptA_LQG_is_initialised(self) result(initialised)
      class(exponential_prop_LQG), intent(in) :: self
      !! Controller
      initialised = self%initialised
   end function exptA_LQG_is_initialised

   subroutine init_exptA_LQG(self, Xdir, Xadj, B, CT, Rinv, Vinv, enable_control)
      class(exponential_prop_LQG), intent(inout)  :: self
      !! Linear Operator.
      class(abstract_sym_low_rank_state_rdp), intent(in) :: Xdir
      !! Riccati matrix representation of solution to direct problem (LQR)
      class(abstract_sym_low_rank_state_rdp), intent(in) :: Xadj
      !! Riccati matrix representation of solution to adjoint problem (LQE)
      class(abstract_vector_rdp), intent(in) :: B(:)
      !! Input matrix
      class(abstract_vector_rdp), intent(in) :: CT(:)
      !! Measurement matrix
      real(dp), intent(in) :: Rinv(:, :)
      !! Inverse control cost
      real(dp), intent(in) :: Vinv(:, :)
      !! Inverse sensor noise
      logical, optional, intent(in) :: enable_control
      !! enable control?
      ! internals
      character(len=*), parameter :: this_procedure = 'init_exptA_LQG'
      ! Initialize propagator.
      call self%prop%initialize(n=2*N, f=rhs_LQG)
      ! setup internal propagator to compute K/L
      call self%prop%setup(Xdir, Xadj, B, CT, Rinv, Vinv)
      ! Choose whether to enable control
      self%control_enabled = optval(enable_control, .true.)
      ! set initialisation flag
      self%initialised = .true.
   end subroutine init_exptA_LQG
   
   subroutine direct_solver_LQG(self, vec_in, vec_out)
      ! Linear Operator.
      class(exponential_prop_LQG), intent(inout)  :: self
      ! Input vector.
      class(abstract_vector_rdp),  intent(in)  :: vec_in
      ! Output vector.
      class(abstract_vector_rdp),  intent(out) :: vec_out
      
      ! Time-integrator.
      character(len=*), parameter :: this_procedure = 'direct_solver_LQG'
      real(dp) :: dt = 1.0_dp
      real(dp) :: state_IC(2*N), state_OC(2*N)
      
      select type(vec_in)
      type is(state_error_vector)
         select type(vec_out)
         type is(state_error_vector)
            ! Get full state + error vector.
            state_IC(:N)   = vec_in%state%state
            state_IC(N+1:) = vec_in%error%state
      
            ! check if control should be enabled
            call self%prop%check_LQG_enabled(self%control_enabled)
            
            ! Integrate forward in time.
            call self%prop%integrate(0.0_dp, state_IC, dt, self%tau, state_OC)
      
            ! Reconstruct state_error_vector.
            vec_out%state%state = state_OC(:N)
            vec_out%error%state = state_OC(N+1:)
         class default
               call type_error('vec_out', 'state_error_vector', 'OUT', this_module, this_procedure)
         end select
      class default
         call type_error('vec_in', 'state_error_vector', 'IN', this_module, this_procedure)
      end select
   end subroutine direct_solver_LQG

end module Ginzburg_Landau_Control