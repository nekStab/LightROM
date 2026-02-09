module Ginzburg_Landau_Tests
   use stdlib_io_npy, only : load_npy
   ! Ginzburg Landau
   use Ginzburg_Landau_Base, only : state_vector
   use Ginzburg_Landau_Operators, only : exponential_prop
   use Ginzburg_Landau_Control, only : exponential_prop_with_control, rks54_class_with_control
   use Ginzburg_Landau_Tests_Lyapunov
   use Ginzburg_Landau_Tests_Riccati
   implicit none

   character(len=*), parameter, private :: this_module = 'Ginzburg_Landau_Tests'

   public :: eigenvalue_analysis
   public :: check_eigenvalues
   public :: integrate_DLRA_fixed
   public :: integrate_DLRA_adaptive

contains

   subroutine eigenvalue_analysis(prop, mold, tmr_name, fname)
      implicit none
      class(abstract_exptA_linop_rdp), intent(inout) :: prop
      class(abstract_vector_rdp), intent(in) :: mold(:)
      character(len=*), optional, intent(in) :: tmr_name
      character(len=*), optional, intent(in) :: fname

      ! internals
      character(len=*), parameter :: this_procedure = 'eigenvalue_analysis'
      character(len=64)                         :: label
      integer                                   :: nev, info
      class(abstract_vector_rdp),           allocatable :: V(:)
      complex(dp),                  allocatable :: lambda(:)
      real(dp),                     allocatable :: residuals(:)

      if (present(tmr_name)) call global_lightROM_timer%add_timer(tmr_name, start=.true.)

      nev = 50; allocate(V(nev), source=mold(1)); call zero_basis(V)
      call eigs(prop, V, lambda, residuals, info, kdim=2*nev)
      call check_info(info, 'eigs', this_module, this_procedure)
      
      lambda = log(lambda)/prop%tau
      if (present(tmr_name)) then
         label = trim(tmr_name)
      else
         label ='eig A'
      end if
      print '(A64,A,F21.7)', padr(trim(label), 64), ': ', maxval(real(lambda))

      if (present(fname)) call save_npy(fname, lambda)
      if (present(tmr_name)) call global_lightROM_timer%stop(tmr_name)
   end subroutine eigenvalue_analysis

   subroutine check_eigenvalues(eq, prop, LTI, U0, Tend, dtv, torder, tolv, adjoint, fname_SVD_base, home)
      implicit none
      character(len=4),              intent(in)    :: eq
      ! propagator
      type(exponential_prop),        intent(inout) :: prop
      ! LTI system
      type(lti_system),              intent(inout) :: LTI
      ! Initial condition
      type(state_vector),            intent(inout) :: U0(:)
      ! End time
      real(dp),                      intent(in)    :: Tend
      ! vector of dt values
      real(dp),                      intent(in)    :: dtv(:)
      ! vector of torders
      integer,                       intent(in)    :: torder
      ! vector of adaptation tolerance
      real(dp),                      intent(in)    :: tolv(:)
      ! number of singular values to print
      logical,                       intent(in)    :: adjoint
      character(len=*),              intent(in)    :: fname_SVD_base
      character(len=128),            intent(in)    :: home

      ! internal
      type(LR_state),                      allocatable :: X
      integer                                          :: j, k
      real(dp)                                         :: tol, tau
      ! Exponential propagator (with control)
      type(rks54_class_with_control),      allocatable :: rkintegrator
      type(exponential_prop_with_control), allocatable :: prop_control
      character(len=256)                               :: fbase, fname, tmr_name, note
      character(len=32)                                :: tolstr, taustr, Tstr
      ! IO
      integer,                              parameter  :: rk_dummy = 10
      logical                                          :: exist_file
      real(dp),                            allocatable :: meta(:)

      ! eig A
      tmr_name =  'eig A'
      fname    = trim(home)//"spectrum_A.npy"
      call eigenvalue_analysis(prop, U0, tmr_name, fname)
      print *, ''

      ! eig A - BK
      X = LR_state()
      if (adjoint) then
         tmr_name = 'eig A-LC: exact'
         fname    = trim(home)//"spectrum_A-LC.npy"
         note = '_adjoint'
      else
         tmr_name = 'eig A-BK: exact'
         fname    = trim(home)//"spectrum_A-BK.npy"
         note = '_direct'
      end if
      call load_X_from_file(X, meta, trim(fname_SVD_base)//eq//trim(note), U0)
      rkintegrator = rks54_class_with_control()
      prop_control = exponential_prop_with_control(1.0_dp, prop=rkintegrator)
      if (adjoint) then
         call prop_control%init(X, LTI%CT, Vinv, adjoint=adjoint, enable_control=.true.)
      else
         call prop_control%init(X, LTI%B, Rinv, adjoint=adjoint, enable_control=.true.)
      end if
      
      call eigenvalue_analysis(prop_control, U0, tmr_name, fname)

      ! deallocate and clean
      deallocate(rkintegrator); deallocate(prop_control)

      print *, ''
      do j = 1, size(tolv)
         tol = tolv(j)
         do k = 1, size(dtv)
            tau = dtv(k)
            note = merge('Padj', 'Pdir', adjoint)
            fbase = make_filename(home, 'DLRA_ADAPT', eq, trim(note), rk_dummy, torder, tau, Tend, tol)
            exist_file = exist_X_file(fbase)
            if (exist_file) then
               ! load X state
               call load_X_from_file(X, meta, fbase, U0)
               ! recreate integrators
               rkintegrator = rks54_class_with_control()
               prop_control = exponential_prop_with_control(1.0_dp, prop=rkintegrator)
               ! initialize integrators
               if (adjoint) then
                  call prop_control%init(X, LTI%CT, Vinv, adjoint=adjoint, enable_control=.true.)
               else
                  call prop_control%init(X, LTI%B, Rinv, adjoint=adjoint, enable_control=.true.)
               end if
               
               call make_labels(Tstr, taustr, tolstr, Tend, tau, tol)
               note = merge('A-LC', 'A-BK', adjoint)
               fname = trim(home)//'spectrum_'//trim(note)//'_Tend'//trim(Tstr)//'_tau'//trim(taustr)//'_tol'//trim(tolstr)//'.npy'
               write(tmr_name,'(*(A))')  'eig ', trim(note), ':   Tend= ', trim(Tstr), '   tau= ', trim(taustr), '   tol= ', trim(tolstr)
               
               ! eigendecomposition
               call eigenvalue_analysis(prop_control, U0, tmr_name, fname)
               
               ! deallocate and clean
               deallocate(rkintegrator); deallocate(prop_control)
            end if
         end do
         write(*,*)
      end do
   end subroutine check_eigenvalues
   
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