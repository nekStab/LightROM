module LightROM_AbstractLTIsystems
   ! Use the abstract linear operator types defined in LightKrylov.
   use LightKrylov, only : abstract_linop_rdp, abstract_vector_rdp, wp => dp
   use LightKrylov_AbstractVectors, only: innerprod
   use LightKrylov_Constants, only : zero_rdp
   implicit none

   private

   character(len=*), parameter :: this_module = 'LR_AbsLTIsys'

   !-------------------------------------------------------
   !-----     ABSTRACT LTI SYSTEM TYPE DEFINITION     -----
   !-------------------------------------------------------

   ! General abstract type for general system.
   type, abstract, public :: abstract_dynamical_system
   end type abstract_dynamical_system

   ! Abstract continuous LTI system.
   type, extends(abstract_dynamical_system), abstract, public :: abstract_lti_system_rdp
      ! Dynamics matrix.
      class(abstract_linop_rdp),  allocatable :: A
      ! Exponential propagator.
      class(abstract_linop_rdp),  allocatable :: prop
      ! Input-to-state matrix.
      class(abstract_vector_rdp), allocatable :: B(:)
      ! State-to-output matrix.
      class(abstract_vector_rdp), allocatable :: CT(:)
      ! Feedthrough matrix.
      real(wp),                   allocatable :: D(:, :)
   contains
   end type abstract_lti_system_rdp

   ! Abstract discrete LTI system.
   type, extends(abstract_dynamical_system), abstract, public :: abstract_dlti_system_rdp
      ! Dynamic matrix.
      class(abstract_linop_rdp),  allocatable :: A
      ! Input-to-state matrix.
      class(abstract_vector_rdp), allocatable :: B(:)
      ! State-to-output matrix.
      class(abstract_vector_rdp), allocatable :: CT(:)
      ! Feedthrough matrix.
      real(wp),                   allocatable :: D(:, :)
      ! Sampling period.
      real(wp)                                :: dt = 1.0_wp
   contains
     private
   end type abstract_dlti_system_rdp

   !--------------------------------------------------------------------
   !-----     ABSTRACT LOW RANK REPRESENTATION TYPE DEFINITION     -----
   !--------------------------------------------------------------------

   ! General abstract type for general system.
   type, abstract, public :: abstract_low_rank_state
   end type abstract_low_rank_state

   type, extends(abstract_low_rank_state), abstract, public :: abstract_low_rank_state_rdp
   end type abstract_low_rank_state_rdp

   ! Abstract symmetric low-rank representation.
   type, extends(abstract_low_rank_state_rdp), abstract, public :: abstract_sym_low_rank_state_rdp
      ! Low-Rank basis.
      class(abstract_vector_rdp),  allocatable :: U(:)
      ! Coefficients
      real(wp),                    allocatable :: S(:, :)
      ! Current approximation rank
      integer                                  :: rk = 1
      ! Simulation time
      real(wp)                                 :: time = 0.0_wp
      ! Total cumulative simulation time
      real(wp)                                 :: tot_time = 0.0_wp
      ! Simulation step counter
      integer                                  :: step = 0
      ! Total cumulative simulation step counter
      integer                                  :: tot_step = 0
      ! Converged?
      logical                                  :: is_converged = .false.
      ! Has rank been initialized? (for rank-adaptive DLRA)
      logical                                  :: rank_is_initialised = .false.
   contains
   end type abstract_sym_low_rank_state_rdp

   !-------------------------------------------------------
   !-----     ABSTRACT CONTROLLER TYPE DEFINITION     -----
   !-------------------------------------------------------

   ! General abstract type for controller.
   type, abstract, public :: abstract_controller
   end type abstract_controller

   type, extends(abstract_controller), abstract, public :: abstract_controller_rdp
   end type abstract_controller_rdp

   ! Abstract Linear-Quadratic Regulator.
   type, extends(abstract_controller_rdp), abstract, public :: abstract_LQR_rdp
      class(abstract_vector_rdp),  allocatable :: U(:)
      !! Low-rank basis for Riccati matrix.
      real(wp),                    allocatable :: RinvBTUS(:,:)
      !! Pre-computed low-rank representation of controller gain.
      integer                                  :: rk_input
      !! Input-to-system rank
      integer                                  :: rk_X
      !! Riccati matrix rank
      logical, private                         :: initialised = .false.
      !! Is initialized
   contains
      ! Initialize controller
      procedure, pass(self), public :: init => init_LQR_rdp
      ! Evaluate control law on input            
      procedure, pass(self), public :: eval => eval_LQR_action_rdp
      ! Check initialisation
      procedure, pass(self), public :: is_initialised => LQR_is_initialised
   end type abstract_LQR_rdp

contains

   !---------------------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE ABSTRACT CONTROLLER TYPES     -----
   !---------------------------------------------------------------------------

   logical function LQR_is_initialised(self) result(is_initialised)
      !! Determine whether the controller has been initialised
      class(abstract_LQR_rdp), intent(in) :: self
      !! Controller
      is_initialised = self%initialised
   end function LQR_is_initialised

   subroutine init_LQR_rdp(self, X, B, Rinv)
      !! Initialise LQR controller and pre-compute action gains
      class(abstract_LQR_rdp),                intent(inout) :: self
      !! Controller
      class(abstract_sym_low_rank_state_rdp), intent(in)    :: X
      !! Low-rank representation of the Riccati matrix
      class(abstract_vector_rdp),             intent(in)    :: B(:)
      !! Input-to-System matrix
      real(wp),                               intent(in)    :: Rinv(:,:)
      !! Inverse control cost.
      ! internals
      real(wp), allocatable :: wrk(:,:)
      ! Sanity check
      call assert_shape(Rinv, [ size(B), size(B) ], 'Rinv', this_module, 'init_LQR_rdp')
      self%rk_input = size(B)
      self%rk_X     = size(X%U)
      allocate(wrk(self%rk_input, self%rk_X), source=zero_rdp)
      allocate(self%U(self%rk_X), source=X%U)
      allocate(self%RinvBTUS(self%rk_input, self%rk_X), source=zero_rdp)
      ! Project input-to-system matrix onto low-rank basis
      wrk = innerprod(B, X%U)                            ! B.T @ U
      ! Multiply with control cost and coefficients
      self%RinvBTUS = -matmul(Rinv, matmul(wrk, X%S))    ! - Rinv @ B.T @ U @ S
      ! set flag
      self%initialised = .true.
   end subroutine init_LQR_rdp

   subroutine eval_LQR_action_rdp(self, action, state)
      !! Evaluate the LQR control law on a state
      class(abstract_LQR_rdp),    intent(in)  :: self
      !! Controller
      real(wp),                   intent(out) :: action(:)
      !! Control action
      class(abstract_vector_rdp), intent(in)  :: state
      !! Current state
      ! internals
      real(wp),                   allocatable :: wrk(:)
      if (.not. self%is_initialised()) then
         call stop_error("The LQR controller is not initialised", this_module, 'eval_LQR_action_rdp')
      end if
      call assert_shape(action, [ self%rk_input ], 'action', this_module, 'eval_LQR_action_rdp')
      ! Project state onto low-rank basis
      allocate(wrk(self%rk_X), source=zero_rdp)
      wrk = innerprod(self%U, state)            ! U.T @ state
      action = matmul(self%RinvBTUS, wrk)       ! - Rinv @ B.T @ U @ S @ U.T @ state
   end subroutine eval_LQR_action_rdp

end module LightROM_AbstractLTIsystems