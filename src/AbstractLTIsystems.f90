module LightROM_AbstractLTIsystems
   use stdlib_optval, only : optval
   ! Use the abstract linear operator types defined in LightKrylov.
   use LightKrylov, only : abstract_linop_rdp, abstract_vector_rdp, dp
   implicit none

   character(len=*), parameter, private :: this_module = 'LR_AbsLTIsys'

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
      real(dp),                   allocatable :: D(:, :)
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
      real(dp),                   allocatable :: D(:, :)
      ! Sampling period.
      real(dp)                                :: dt = 1.0_dp
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
      real(dp),                    allocatable :: S(:, :)
      ! Current approximation rank
      integer                                  :: rk = 1
      ! Simulation time
      real(dp)                                 :: time = 0.0_dp
      ! Total cumulative simulation time
      real(dp)                                 :: tot_time = 0.0_dp
      ! Simulation step counter
      integer                                  :: step = 0
      ! Total cumulative simulation step counter
      integer                                  :: tot_step = 0
      ! Converged?
      logical                                  :: is_converged = .false.
      ! Has rank been initialized? (for rank-adaptive DLRA)
      logical                                  :: rank_is_initialised = .false.
   contains
      procedure, pass(self), public :: reset => abstract_sym_low_rank_state_reset
      procedure, pass(self), public :: increment_counters => abstract_sym_low_rank_state_increment_counters
   end type abstract_sym_low_rank_state_rdp

contains

   subroutine abstract_sym_low_rank_state_reset(self, full, tot_time, tot_step)
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: self
      logical, optional, intent(in) :: full
      real(dp), optional, intent(in) :: tot_time
      integer, optional, intent(in) :: tot_step
      ! internal
      logical :: full_
      full_ = optval(full, .false.)
      self%time = 0.0_dp
      self%step = 0
      self%is_converged = .false.
      if (full_) then
         self%rank_is_initialised = .false.
         self%tot_time = 0.0_dp
         self%tot_step = 0
      end if
      if (present(tot_time)) self%tot_time = tot_time
      if (present(tot_step)) self%tot_step = tot_step
   end subroutine abstract_sym_low_rank_state_reset

   subroutine abstract_sym_low_rank_state_increment_counters(self, tau)
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: self
      real(dp), intent(in) :: tau
      self%time = self%time + tau
      self%step = self%step + 1
      self%tot_time = self%tot_time + tau
      self%tot_step = self%tot_step + 1
   end subroutine abstract_sym_low_rank_state_increment_counters
end module LightROM_AbstractLTIsystems