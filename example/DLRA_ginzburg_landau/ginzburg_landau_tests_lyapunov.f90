module Ginzburg_Landau_Tests_Lyapunov
   ! Standard Library.
   use stdlib_optval, only : optval
   use stdlib_linalg, only : diag, svd, svdvals
   use stdlib_io_npy, only : save_npy, load_npy
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov_AbstractVectors ! linear_combination
   use LightKrylov_Utils ! svd, sqrtm
   ! LightROM
   use LightROM_Utils ! Balancing_Transformation
   ! Lyapunov Solver
   use LightROM_LyapunovSolvers
   use LightROM_LyapunovUtils
   ! Ginzburg Landau
   use Ginzburg_Landau_Base
   use Ginzburg_Landau_Operators
   use Ginzburg_Landau_RKlib
   use Ginzburg_Landau_Utils
   implicit none

   character(len=*), parameter, private :: this_module = 'Ginzburg_Landau_Tests_Lyapunov'
   ! reference solutiions using RK
   character(len=128), parameter :: fname_RK_base = 'Xref_RK_'

contains

   !-------------------------------------------------------------------------------------------
   !
   ! LYAPUNOV EQUATION
   !
   !-------------------------------------------------------------------------------------------

   subroutine run_lyapunov_reference_RK(LTI, Xref, Xref_RK, U0, S0, Tend, nrep, iref, adjoint, home, if_save_output)
      ! LTI system
      type(lti_system),              intent(inout) :: LTI
      ! Reference solution (BS)
      real(dp),                      intent(in)    :: Xref(N,N)
      ! Reference solution (RK)
      real(dp),                      intent(out)   :: Xref_RK(N,N)
      ! Initial condition
      type(state_vector),            intent(inout) :: U0(:)
      real(dp),                      intent(inout) :: S0(:,:)
      real(dp),                      intent(in)    :: Tend
      integer,                       intent(in)    :: nrep
      integer,                       intent(in)    :: iref
      logical,                       intent(in)    :: adjoint
      character(len=128),            intent(in)    :: home
      logical,                       intent(in)    :: if_save_output

      ! Internals
      type(rk_lyapunov),             allocatable   :: RK_propagator
      type(state_matrix)                           :: X_mat(2)
      real(dp),                      allocatable   :: X_RK(:,:,:), U_load(:,:), meta(:)
      integer                                      :: irep
      real(dp)                                     :: etime0, etime, Tstep
      ! OUTPUT
      real(dp)                                     :: X_out(N,N)
      character(len=256)                           :: note, fbase
      logical                                      :: exist_file
      ! timing
      integer                                      :: clock_rate, clock_start, clock_stop

      character(len=*), parameter :: case = 'RKLIB'
      character(len=*), parameter :: eq   = 'lyap'

      call system_clock(count_rate=clock_rate)

      print *, ''
      print *, '#######################################################'
      print *, '#                                                     #'
      print *, '#             Solution using Runge-Kutta              #'
      print *, '#                                                     #'
      print *, '#######################################################'
      print *, ''
      print '(2X,A,F6.2)', 'Tend = ', Tend
      print '(2X,A,I3)', 'nrep = ', nrep
      print '(2X,A,I3)', 'iref = ', iref
      print '(2X,A,L3)', 'save = ', if_save_output
      print *, ''

      call print_header(case, eq)
      write(*,'(I7,F10.4,3(1X,E18.6),F10.4," s",A)') 0, 0.0, norm2(X_out)/N, norm2(X_out - Xref)/N, &
                                 & norm2(CALE(X_out, adjoint))/N, 0.0, ''

      ! init
      allocate(X_RK(N, N, nrep))
      Tstep = Tend/nrep
      ! initialize exponential propagator
      RK_propagator = RK_lyapunov(Tstep)
      ! Set initial condition for RK
      call reconstruct_solution(X_out, U0, S0)
      call set_state(X_mat(1:1), X_out, 'Set initial condition')

      fbase = make_filename_RK(home, fname_RK_base, eq, Tend, adjoint)
      exist_file = exist_RK_file(fbase)
      if (exist_file) then
         call load_npy(trim(fbase)//'_X.npy', U_load)
         call load_npy(trim(fbase)//'_meta.npy', meta)
         irep  = 1
         X_RK(:,:,irep) = U_load
         etime = meta(1)
         note = '   < reference'
         call print_rklib_output(eq, irep, Tend, X_RK, Xref, etime, note, adjoint)
         Xref_RK(:,:) = X_RK(:,:,irep)
      else
         etime = 0.0_dp
         do irep = 1, nrep
            call system_clock(count=clock_start)     ! Start Timer
            ! integrate
            if (adjoint) then
               call RK_propagator%rmatvec(X_mat(1), X_mat(2))
            else
               call RK_propagator%matvec(X_mat(1), X_mat(2))
            end if
            call system_clock(count=clock_stop)      ! Stop Timer
            etime0 = real(clock_stop-clock_start)/real(clock_rate)
            etime = etime + etime0
            ! recover output
            call get_state(X_RK(:,:,irep), X_mat(2:2), 'Extract RK solution')
            ! replace input
            call set_state(X_mat(1:1), X_RK(:,:,irep), 'Reset initial condition')
            ! print information
            note = merge('   < reference', '              ', irep == iref)
            call print_rklib_output(eq, irep, Tstep, X_RK, Xref, etime0, note, adjoint)
         enddo
         Xref_RK(:,:) = X_RK(:,:,iref)
         call get_metadata(meta, eq, 0, 0, 0.0_dp, 0, Tend, Xref_RK, Xref_RK, Xref_RK, etime, adjoint)
      end if
      ! save output
      if (if_save_output .and. .not. exist_file) call save_RK_state_npy(fbase, Xref_RK, meta)
      print *, ''
      print '(A,F16.12)', '  |  X_RK  |/N = ', norm2(Xref_RK)/N
      print '(A,F16.12)', '  | res_RK |/N = ', norm2(CALE(Xref_RK, adjoint))/N
      return
   end subroutine run_lyapunov_reference_RK

   subroutine run_lyapunov_DLRA_test(LTI, Xref, Xref_RK, U0, S0, Tend, dtv, rkv, TOv, nprint, adjoint, home, if_save_output)
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
      ! Adjoint mode?
      logical,                       intent(in)    :: adjoint
      ! Home folder
      character(len=128),            intent(in)    :: home
      ! save data to python format?
      logical,                       intent(in)    :: if_save_output

      ! Internals
      type(LR_state),                allocatable   :: X
      integer                                      :: info, i, j, k, rk, torder, irep, nrep, nsteps
      real(dp)                                     :: etime, tau, Tstep, Ttot
      ! Settings
      type(dlra_opts)                              :: opts
      logical                                      :: expta_adj
      ! IO
      character(len=256)                           :: fbase, fchomp
      real(dp)                                     :: X_out(N,N)
      logical                                      :: exist_file
      real(wp),                      allocatable   :: meta(:)
      character(len=128)                           :: note
      ! timing
      integer                                      :: clock_rate, clock_start, clock_stop
      
      character(len=*), parameter :: case = 'DLRA_FIXED'
      character(len=*), parameter :: eq   = 'lyap'
      
      call system_clock(count_rate=clock_rate)

      print *, ''
      print *, '#######################################################'
      print *, '#                                                     #'
      print *, '#           Solution using fixed-rank DLRA            #'
      print *, '#                                                     #'
      print *, '#######################################################'
      print *, ''
      print '(2X,A6,*(I6,1X))',    padr('TO  =',6), TOv
      print '(2X,A6,*(I6,1X))',    padr('rk  =',6), rkv
      print '(2X,A6,*(ES6.0,1X))', padr('tau =',6), dtv
      print *, ''
      call print_header(case, eq)

      ! init
      X = LR_state()
      Tstep = one_rdp
      note = merge('Yobs', 'Xctl', adjoint)
      expta_adj = adjoint     ! the adjoint Lyapunov equation uses the adjoint expta solver
      opts = dlra_opts(chktime=one_rdp, inc_tol=atol_dp, if_rank_adaptive=.false.)

      do i = 1, size(rkv)
         rk = rkv(i)
         do j = 1, size(TOv)
            torder = TOv(j)
            do k = 1, size(dtv)
               tau = dtv(k)
               nsteps = nint(Tend/tau)
               
               ! set solver options
               opts%mode = torder

               fbase = make_filename(home, case, eq, note, rk, torder, tau, Tend)
               fchomp = replace_all(fbase, trim(home), '')
               exist_file = exist_X_file(fbase)
               if (exist_file) then
                  call load_X_from_file(X, meta, fbase, U0)
                  etime = meta(1)
               else
                  call reset_timers()
                  ! Initialize low-rank representation with rank rk
                  call X%initialize_LR_state(U0, S0, rk, rkmax, .false.)
                  ! run integrator
                  call system_clock(count=clock_start)     ! Start Timer
                  if (adjoint) then
                     call projector_splitting_DLRA_lyapunov_integrator(X, LTI%prop, LTI%CT, Tend, tau, info, &
                                                               & exptA=exptA, iftrans=expta_adj, options=opts)
                  else
                     call projector_splitting_DLRA_lyapunov_integrator(X, LTI%prop, LTI%B, Tend, tau, info, &
                                                               & exptA=exptA, iftrans=expta_adj, options=opts)
                  end if
                  call system_clock(count=clock_stop)      ! Stop Timer
                  call reset_lyapunov_solver(home, fchomp)
                  etime = real(clock_stop-clock_start)/real(clock_rate)
               end if
               ! Reconstruct solution
               call reconstruct_solution(X_out, X)
               ! Compute metadata
               call get_metadata(meta, eq, rk, torder, tau, nsteps, tend, X_out, Xref_RK, Xref, etime, adjoint)
               ! print information
               call print_dlra_output(eq, note, meta)
               ! save output
               if (if_save_output .and. .not. exist_file) then
                  call save_LR_state_npy(fbase, X, weight_mat, meta)
                  call save_metadata(fbase, case, eq, rk, torder, tau, nsteps, Tend, adjoint)
               end if
               ! cleanup
               deallocate(X%U); deallocate(X%S)
            end do
            print *, ''
         end do
      end do
      if (nprint > 0) then
         call print_svdvals(Xref,    'X_BS', nprint)
         call print_svdvals(Xref_RK, 'X_RK', nprint)
         call print_svdvals(X_out,   'X_D ', nprint)
      end if

   end subroutine run_lyapunov_DLRA_test

   subroutine run_lyapunov_DLRArk_test(LTI, Xref, Xref_RK, U0, S0, Tend, dtv, TOv, tolv, nprint, adjoint, home, if_save_output)
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
      ! Adjoint mode?
      logical,                       intent(in)    :: adjoint
      ! Home folder
      character(len=128),            intent(in)    :: home
      ! save data to python format?
      logical,                       intent(in)    :: if_save_output
      
      ! Internals
      type(LR_state),                allocatable   :: X
      integer                                      :: info, i, j, k, rk, torder, nsteps
      real(dp)                                     :: etime, tau
      ! Settings
      type(dlra_opts)                              :: opts
      logical                                      :: expta_adj
      ! IO
      character(len=256)                           :: fbase, fchomp
      real(dp)                                     :: X_out(N,N)
      logical                                      :: exist_file
      real(wp),                      allocatable   :: meta(:)
      character(len=128)                           :: note
      ! timing
      integer                                      :: clock_rate, clock_start, clock_stop

      character(len=*), parameter :: case = 'DLRA_ADAPT'
      character(len=*), parameter :: eq   = 'lyap'
      
      call system_clock(count_rate=clock_rate)
      
      print *, ''
      print *, '#######################################################'
      print *, '#                                                     #'
      print *, '#         Solution using rank-adaptive DLRA           #'
      print *, '#                                                     #'
      print *, '#######################################################'
      print *, ''
      print '(2X,A6,*(I6,1X))',    padr('TO  =',6), TOv
      print '(2X,A6,*(ES6.0,1X))', padr('tol =',6), tolv
      print '(2X,A6,*(ES6.0,1X))', padr('tau =',6), dtv
      print *, ''
 
      ! init
      X = LR_state()
      rk = rk_X0_lyapunov
      note = merge('Yobs', 'Xctl', adjoint)
      expta_adj = adjoint     ! the adjoint Lyapunov equation uses the adjoint expta solver
      opts = dlra_opts(chktime=one_rdp, inc_tol=atol_dp, if_rank_adaptive=.true.)

      do i = 1, size(tolv)
         call print_header(case, eq)
         opts%tol = tolv(i)
         print '(A,E9.2)', ' SVD tol = ', opts%tol
         print *, ''
         do j = 1, size(TOv)
            torder = TOv(j)
            do k = 1, size(dtv)
               tau = dtv(k)
               nsteps = nint(Tend/tau)
               
               ! set solver options
               opts%mode = torder

               fbase = make_filename(home, case, eq, note, rk, torder, tau, Tend, opts%tol)
               fchomp = replace_all(fbase, trim(home), '')
               exist_file = exist_X_file(fbase)
               if (exist_file) then
                  call load_X_from_file(X, meta, fbase, U0)
                  etime = meta(1)
               else
                  call reset_timers()
                  ! Initialize low-rank representation with rank rk
                  call X%initialize_LR_state(U0, S0, rk, rkmax, opts%if_rank_adaptive)
                  ! run integrator
                  call system_clock(count=clock_start)     ! Start Timer
                  if (adjoint) then
                     call projector_splitting_DLRA_lyapunov_integrator(X, LTI%prop, LTI%CT, Tend, tau, info, &
                                                                  & exptA=exptA, iftrans=.true., options=opts)
                  else
                     call projector_splitting_DLRA_lyapunov_integrator(X, LTI%prop, LTI%B, Tend, tau, info, &
                                                                  & exptA=exptA, iftrans=.false., options=opts)
                  end if
                  call system_clock(count=clock_stop)      ! Stop Timer
                  call reset_lyapunov_solver(home, fchomp)
                  etime = real(clock_stop-clock_start)/real(clock_rate)
               end if
               ! Reconstruct solution
               call reconstruct_solution(X_out, X)
               ! Compute metadata
               call get_metadata(meta, eq, X%rk, torder, tau, nsteps, tend, X_out, Xref_RK, Xref, etime, adjoint)
               ! print information
               call print_dlra_output(eq, note, meta)
               ! save output
               if (if_save_output .and. .not. exist_file) then
                  call save_LR_state_npy(fbase, X, weight_mat, meta)
                  call save_metadata(fbase, case, eq, X%rk, torder, tau, nsteps, Tend, adjoint)
               end if
               ! cleanup
               deallocate(X%U); deallocate(X%S)
            end do
            print *, ''
         end do
         call print_svdvals(Xref,    'X_BS', nprint)
         call print_svdvals(Xref_RK, 'X_RK', nprint)
         call print_svdvals(X_out,   'X_D ', nprint)
      end do

   end subroutine run_lyapunov_DLRArk_test

end module Ginzburg_Landau_Tests_Lyapunov