module LightROM_Timing
   use LightKrylov
   use LightKrylov, only : dp, wp => dp
   use LightKrylov_Logger
   use LightKrylov_Timing
   implicit none
   private
   character(len=*), parameter :: this_module      = 'LR_Timer'
   character(len=*), parameter :: this_module_long = 'LightROM_Timer'
   logical :: if_time = .false.

   public :: time_lightrom
   public :: global_lightROM_timer

   ! LightROM_watch type
   type, extends(abstract_watch), public :: lightROM_watch
      !! Global timing structure to contain all timers within LightROM
      private
      character(len=128) :: name = 'lightROM_timer'
   contains
      private
      procedure, pass(self), public :: set_private_timers => set_lightROM_timers
   end type lightROM_watch

   type(lightROM_watch) :: global_lightROM_timer

contains

   logical function time_lightROM() result(if_time_lightROM)
      if_time_lightROM = if_time
   end function time_lightROM

   !--------------------------------------------------------------
   !  Concrete implementations for the lightROM_watch type
   !--------------------------------------------------------------

   subroutine set_lightROM_timers(self)
      !! Initialize global watch within LightROM and define private system timers.
      class(lightROM_watch), intent(inout) :: self
      ! internal
      integer :: istart
      ! timers for LightKrylov_LyapunovSolvers
      istart = self%get_timer_count() + 1
      call self%add_timer('projector_splitting_DLRA_lyapunov_integrator_rdp')
      call self%add_timer('projector_splitting_DLRA_lyapunov_step_rdp')
      call self%add_timer('M_forward_map_rdp')
      call self%add_timer('G_forward_map_rdp')
      call self%add_timer('K_step_lyapunov_rdp')
      call self%add_timer('S_step_lyapunov_rdp')
      call self%add_timer('L_step_lyapunov_rdp')
      call self%add_group('LyapunovSolvers', istart=istart, iend=self%get_timer_count())
      ! timers for LightKrylov_RiccatiSolvers
      istart = self%get_timer_count() + 1
      call self%add_timer('projector_splitting_DLRA_riccati_integrator_rdp')
      call self%add_timer('projector_splitting_DLRA_riccati_step_rdp')
      call self%add_timer('G_forward_map_riccati_rdp')
      call self%add_timer('K_step_lyapunov_riccati_rdp')
      call self%add_timer('S_step_lyapunov_riccati_rdp')
      call self%add_timer('L_step_lyapunov_riccati_rdp')
      call self%add_group('RiccatiSolvers', istart=istart, iend=self%get_timer_count())
   end subroutine set_lightROM_timers

end module LightROM_Timing
