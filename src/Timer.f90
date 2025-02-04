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

   public :: initialize_timers
   public :: reset_timers
   public :: enumerate_timers
   public :: finalize_timers

   ! LightROM_watch type
   type, extends(abstract_watch), public :: lightROM_watch
      !! Global timing structure to contain all timers within LightROM
   contains
      private
      procedure, pass(self), public :: set_private_timers_and_name => set_lightROM_timers
   end type lightROM_watch

   type(lightROM_watch) :: global_lightROM_timer

contains

   logical function time_lightROM() result(if_time_lightROM)
      if_time_lightROM = if_time
   end function time_lightROM

   subroutine set_lightROM_timer_switch(value)
      logical, intent(in) :: value     
      if (if_time .neqv. value) then
         if_time = value
         if (if_time) then
            call log_message('LightROM timing enabled.', module=this_module)
         else
            call log_message('LightROM timing disabled.', module=this_module)
         end if
      else
         call log_debug('LightROM timing switched unchanged.', module=this_module)
      end if      
   end subroutine set_lightROM_timer_switch

   !--------------------------------------------------------------
   !  Concrete implementations for the lightROM_watch type
   !--------------------------------------------------------------

   subroutine set_lightROM_timers(self)
      !! Initialize global watch within LightROM and define private system timers.
      class(lightROM_watch), intent(inout) :: self
      ! internal
      integer :: istart, iend
      call self%set_watch_name('LightROM_timer')
      ! timers for LightKrylov_LyapunovSolvers
      call self%add_timer('DLRA_lyapunov_integrator_rdp', count=istart)
      call self%add_timer('DLRA_lyapunov_step_rdp')
      call self%add_timer('M_forward_map_rdp')
      call self%add_timer('G_forward_map_lyapunov_rdp')
      call self%add_timer('K_step_lyapunov_rdp')
      call self%add_timer('S_step_lyapunov_rdp')
      call self%add_timer('L_step_lyapunov_rdp', count=iend)
      ! define LyapunovSolvers group
      call self%add_group('LyapunovSolvers', istart=istart, iend=iend)
      ! timers for LightKrylov_RiccatiSolvers
      call self%add_timer('DLRA_riccati_integrator_rdp', count=istart)
      call self%add_timer('DLRA_riccati_step_rdp')
      call self%add_timer('G_forward_map_riccati_rdp')
      call self%add_timer('K_step_riccati_rdp')
      call self%add_timer('S_step_riccati_rdp')
      call self%add_timer('L_step_riccati_rdp', count=iend)
      ! define RiccatiSolvers group
      call self%add_group('RiccatiSolvers', istart=istart, iend=iend)
      ! Enable timing
      call set_lightROM_timer_switch(.true.)
   end subroutine set_lightROM_timers

   subroutine initialize_timers
      call global_lightkrylov_timer%initialize()
      call global_lightROM_timer%initialize()
   end subroutine initialize_timers

   subroutine reset_timers(soft, clean)
      logical, optional, intent(in) :: soft
      logical, optional, intent(in) :: clean
      call global_lightkrylov_timer%reset_all(soft, clean)
      call global_lightROM_timer%reset_all(soft, clean)
   end subroutine reset_timers

   subroutine enumerate_timers(only_user)
      logical, optional, intent(in) :: only_user
      call global_lightkrylov_timer%enumerate(only_user)
      call global_lightROM_timer%enumerate(only_user)
   end subroutine enumerate_timers

   subroutine finalize_timers
      call global_lightkrylov_timer%finalize()
      call global_lightROM_timer%finalize()
   end subroutine finalize_timers

end module LightROM_Timing
