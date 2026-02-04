module Ginzburg_Landau_Tests
   use stdlib_io_npy, only : load_npy
   ! Ginzburg Landau
   use Ginzburg_Landau_Base, only : state_vector
   use Ginzburg_Landau_Operators, only : exponential_prop
   use Ginzburg_Landau_Control, only : exponential_prop_with_control
   use Ginzburg_Landau_Tests_Lyapunov
   use Ginzburg_Landau_Tests_Riccati
   implicit none

   character(len=*), parameter, private :: this_module = 'Ginzburg_Landau_Tests'

   public :: eigenvalue_analysis
   public :: eigenvalue_analysis_control
   public :: integrate_DLRA_fixed
   public :: integrate_DLRA_adaptive

contains

   subroutine eigenvalue_analysis(prop, tmr_name, fname)
      implicit none
      type(exponential_prop), intent(inout) :: prop
      character(len=*), optional, intent(in) :: tmr_name
      character(len=*), optional, intent(in) :: fname

      ! internals
      character(len=*), parameter :: this_procedure = 'eigenvalue_analysis'
      integer                                   :: nev, info
      type(state_vector),           allocatable :: V(:)
      complex(dp),                  allocatable :: lambda(:)
      real(dp),                     allocatable :: residuals(:)

      if (present(tmr_name)) call global_lightROM_timer%add_timer(tmr_name, start=.true.)

      nev = 50; allocate(V(nev)); call zero_basis(V)
      call eigs(prop, V, lambda, residuals, info, kdim=2*nev)
      call check_info(info, 'eigs', this_module, this_procedure)
      
      lambda = log(lambda)/prop%tau
      print '(A64,A,F21.7)', padr('eig A', 64), ': ', maxval(real(lambda))

      if (present(fname)) call save_npy(fname, lambda)
      if (present(tmr_name)) call global_lightROM_timer%stop(tmr_name)
   end subroutine eigenvalue_analysis

   subroutine eigenvalue_analysis_control(prop, tmr_name, fname)
      implicit none
      type(exponential_prop_with_control), intent(inout) :: prop
      character(len=*), optional, intent(in) :: tmr_name
      character(len=*), optional, intent(in) :: fname

      ! internals
      character(len=*), parameter :: this_procedure = 'eigenvalue_analysis_control'
      integer                                   :: nev, info
      type(state_vector),           allocatable :: V(:)
      complex(dp),                  allocatable :: lambda(:)
      real(dp),                     allocatable :: residuals(:)

      if (present(tmr_name)) call global_lightROM_timer%add_timer(tmr_name, start=.true.)

      nev = 50; allocate(V(nev)); call zero_basis(V)
      call eigs(prop, V, lambda, residuals, info, kdim=2*nev)
      call check_info(info, 'eigs', this_module, this_procedure)
      
      lambda = log(lambda)/prop%tau
      print '(A64,A,F21.7)', padr('eig A-BK', 64), ': ', maxval(real(lambda))

      if (present(fname)) call save_npy(fname, lambda)
      if (present(tmr_name)) call global_lightROM_timer%stop(tmr_name)
   end subroutine eigenvalue_analysis_control
   
   subroutine integrate_RK(eq, LTI, Xref, Xref_RK, U0, S0, Tend, adjoint, home, main_run)
      implicit none
      character(len=4),              intent(in)    :: eq
      ! LTI system
      type(lti_system),              intent(inout) :: LTI
      ! Reference solution (SD)
      real(dp),                      intent(in)    :: Xref(N,N)
      ! Reference solution (RK)
      real(dp),                      intent(out)   :: Xref_RK(N,N)
      ! Initial condition
      type(state_vector),            intent(inout) :: U0(:)
      real(dp),                      intent(inout) :: S0(:,:)
      real(dp),                      intent(in)    :: Tend
      logical,                       intent(in)    :: adjoint
      character(len=128),            intent(in)    :: home
      logical,                       intent(in)    :: main_run

      ! internal
      character(len=128) :: tmr_name
      logical :: if_save_output
      integer :: nstep, iref
      
      if (main_run) then
         tmr_name = 'Runge-Kutta'
         if_save_output = .true.
         nstep = int(Tend)
         iref  = nstep
      else
         tmr_name = 'Short time: Runge-Kutta'
         if_save_output = .false.
         nstep = int(Tend/0.001_dp)
         iref  = nstep
      end if
      
      call global_lightROM_timer%add_timer(tmr_name, start=.true.)
      if (eq == 'lyap') then
         call run_lyapunov_reference_RK(LTI, Xref, Xref_RK, U0, S0, Tend, nstep, iref, adjoint, home, if_save_output)
      else
         call run_riccati_reference_RK( LTI, Xref, Xref_RK, U0, S0, Tend, nstep, iref, adjoint, home, if_save_output)
      end if
      call global_lightROM_timer%stop(tmr_name)
   end subroutine integrate_RK
   
   subroutine integrate_DLRA_fixed(eq, LTI, Xref, Xref_RK, U0, S0, Tend, dtv, rkv, TOv, nprint, adjoint, home, main_run)
      implicit none
      character(len=4),              intent(in)    :: eq
      ! LTI system
      type(lti_system),              intent(inout) :: LTI
      ! Reference solution (BS)
      real(dp),                      intent(in)    :: Xref(N,N)
      ! Reference solution (RK)
      real(dp),                      intent(in)    :: Xref_RK(N,N)
      ! Initial condition
      type(state_vector),            intent(inout) :: U0(:)
      real(dp),                      intent(inout) :: S0(:,:)
      real(dp),                      intent(in)    :: Tend
      ! vector of dt values
      real(dp),                      intent(in)    :: dtv(:)
      ! vector of rank values
      integer,                       intent(in)    :: rkv(:)
      ! vector of torders
      integer,                       intent(in)    :: TOv(:)
      ! number of singular values to print
      integer,                       intent(in)    :: nprint
      logical,                       intent(in)    :: adjoint
      character(len=128),            intent(in)    :: home
      logical,                       intent(in)    :: main_run
      
      ! internal
      character(len=128) :: tmr_name
      logical :: if_save_output

      if (main_run) then
         if_save_output = .true. 
         tmr_name = 'fixed-rank DLRA'
      else
         if_save_output = .false. 
         tmr_name = 'Short_time: fixed-rank DLRA'
      end if
      call global_lightROM_timer%add_timer(tmr_name, start=.true.)
      if (eq == 'lyap') then
         call run_lyapunov_DLRA_test(LTI, Xref, Xref_RK, U0, S0, Tend, dtv, rkv, TOv, nprint, adjoint, home, if_save_output)
      else         
         call run_riccati_DLRA_test(LTI, Xref, Xref_RK, U0, S0, Tend, dtv, rkv, TOv, nprint, adjoint, home, if_save_output)
      end if
      call global_lightROM_timer%stop(tmr_name)
   end subroutine integrate_DLRA_fixed

   subroutine integrate_DLRA_adaptive(eq, LTI, Xref, Xref_RK, U0, S0, Tend, dtv, TOv, tolv, nprint, adjoint, home, main_run)
      implicit none
      character(len=4),              intent(in)    :: eq
      ! LTI system
      type(lti_system),              intent(inout) :: LTI
      ! Reference solution (BS)
      real(dp),                      intent(in)    :: Xref(N,N)
      ! Reference solution (RK)
      real(dp),                      intent(in)    :: Xref_RK(N,N)
      ! Initial condition
      type(state_vector),            intent(inout) :: U0(:)
      real(dp),                      intent(inout) :: S0(:,:)
      real(dp),                      intent(in)    :: Tend
      ! vector of dt values
      real(dp),                      intent(in)    :: dtv(:)
      ! vector of torders
      integer,                       intent(in)    :: TOv(:)
      ! vector of adaptation tolerance
      real(dp),                      intent(in)    :: tolv(:)
      ! number of singular values to print
      integer,                       intent(in)    :: nprint
      logical,                       intent(in)    :: adjoint
      character(len=128),            intent(in)    :: home
      logical,                       intent(in)    :: main_run

      ! internal
      character(len=128) :: tmr_name
      logical :: if_save_output

      if (main_run) then
         if_save_output = .true. 
         tmr_name = 'rank-adaptive DLRA'
      else
         if_save_output = .false. 
         tmr_name = 'Short_time: rank-adaptive DLRA'
      end if
      call global_lightROM_timer%add_timer(tmr_name, start=.true.)
      if (eq == 'lyap') then
         call run_lyapunov_DLRArk_test(LTI, Xref, Xref_RK, U0, S0, Tend, dtv, TOv, tolv, nprint, adjoint, home, if_save_output)
      else
         call run_riccati_DLRArk_test(LTI, Xref, Xref_RK, U0, S0, Tend, dtv, TOv, tolv, nprint, adjoint, home, if_save_output)
      end if
      call global_lightROM_timer%stop(tmr_name)

   end subroutine integrate_DLRA_adaptive

end module Ginzburg_Landau_Tests