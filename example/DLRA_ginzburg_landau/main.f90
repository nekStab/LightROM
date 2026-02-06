program demo
   ! Standard Library.
   use stdlib_strings, only: padr
   use stdlib_optval, only : optval 
   use stdlib_linalg, only : eye, diag
   use stdlib_math, only : all_close, logspace
   use stdlib_io_npy, only : save_npy, load_npy
   use stdlib_logger, only : information_level, warning_level, debug_level, error_level, none_level
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov_Logger
   use LightKrylov_AbstractVectors
   use LightKrylov_ExpmLib
   use LightKrylov_Utils
   ! LightROM
   use LightROM_AbstractLTIsystems
   use LightROM_Utils
   use LightROM_Timing
   use LightROM_LyapunovSolvers
   use LightROM_LyapunovUtils
   ! GInzburg-Landau
   use Ginzburg_Landau_Base
   use Ginzburg_Landau_Operators
   use Ginzburg_Landau_Control
   use Ginzburg_Landau_Utils
   use Ginzburg_Landau_Tests_Lyapunov
   use Ginzburg_Landau_Tests_Riccati
   use Ginzburg_Landau_Tests
   implicit none

   character(len=*), parameter :: this_module = 'Ginzburg_Landau_Main'

   real(dp), parameter :: Tend_long  = 50.0
   real(dp), parameter :: Tend_short = 0.01
   character(len=128), parameter :: home_base = 'example/DLRA_ginzburg_landau/local/'
   ! reference solutiions using SVD
   character(len=128), parameter :: fname_SVD_base = trim(home_base)//'Xref_SVD_'
   ! rk_B & rk_C are set in ginzburg_landau_base.f90

   integer  :: nrk, ntau, rk,  torder
   real(dp) :: tau, T_POD, Tend
   ! vector of dt values
   real(dp), allocatable :: dtv(:), tolv(:), Tendv(:)
   ! vector of rank values
   integer, allocatable :: rkv(:), TOv(:)

   ! Exponential propagator (RKlib).
   type(GL_operator),            allocatable :: A
   type(exponential_prop),       allocatable :: prop

   ! Exponential propagator (with control)
   type(rks54_class_with_control), allocatable :: rkintegrator
   type(exponential_prop_with_control), allocatable :: prop_control

   ! LTI system
   type(lti_system)                          :: LTI
   type(dlra_opts)                           :: opts
   
   type(LR_state),               allocatable :: X

   ! Initial condition
   type(state_vector),           allocatable :: U0(:)!, output(:)
   real(dp),                     allocatable :: S0(:,:)
   
   ! OUTPUT
   real(dp)                                  :: X_out(N,N)

   ! Reference solutions (BS & RK)
   real(dp)                                  :: Xref(N,N)
   real(dp)                                  :: Xref_RK(N,N)

   ! IO
   real(dp),                     allocatable :: meta(:)
   
   ! POD
   type(state_vector),           allocatable :: X0(:)

   ! Information flag.
   integer                                   :: info

   ! Misc
   integer                                   :: i, j, k, is, ie, it, irep
   integer                                   :: rk_X0
   character(len=2)                          :: refid
   character(len=4)                          :: eq
   character(len=32)                         :: tolstr, taustr, Tstr
   character(len=256)                        :: fbase, fname, tmr_name, note
   real(dp),                     allocatable :: svals(:)
   real(dp)                                  :: tol
   integer                                   :: nprint
   logical                                   :: exist_file
   character(len=128)                        :: home

   !--------------------------------
   ! Define which examples to run:
   !
   logical, parameter :: if_lyapunov = .false.
   !
   ! if_lyapunov = .true.:  Solve the Lyapunov equation:   0 = A @ X + X @ A.T + Q
   !
   ! if_lyapunov = .false.: Solve the Riccati equation:    0 = A @ X + X @ A.T + X @ B @ @ R^{-1} @ B.T @ W @ X + Q
   !
   logical, parameter :: if_adj = .false.
   ! Only considered if if_lyapunov = .true.
   !
   ! Adjoint = .true.:      Solve the adjoint Lyapunov equation:  0 = A.T @ X + X @ A + C.T @ C @ W
   !     The solution to this equation is called the observability Gramian Y.
   !
   ! Adjoint = .false.:     Solve the direct Lyapunov equation:   0 = A @ X + X @ A.T + B @ B.T @ W
   !     The solution to this equation is called the controllability Gramian X.
   !
   logical, parameter :: main_run = .true.
   !
   ! Run the computation instead of the test
   !
   logical, parameter :: run_fixed_rank_test = .true.
   !
   ! Integrate the same initial condition to steady state with Runge-Kutta and DLRA.
   !
   ! As the steady state is approached, the error/residual for Runge-Kutta goes to zero.
   ! Similarly, the test shows the effect of step size, rank and temporal order on the solution
   ! using DLRA.
   !
   logical, parameter :: run_rank_adaptive_test = .true.
   !
   ! Integrate the same initial condition to steady state with Runge-Kutta and DLRA using an 
   ! adaptive rank.
   !
   ! The DLRA algorthm automatically determines the rank necessary to integrate the equations
   ! such that the error on the singular values does not exceed a chosen tolerance. This rank
   ! depends on the tolerance but also the chosen time-step.
   !
   logical, parameter :: run_eigenvalue_test = .true.
   !
   ! Check the control efficacy using the eigenvalues of the controlled system matrix
   !
   !--------------------------------
   Tend = merge(Tend_long, Tend_short, main_run)
   eq   = merge('lyap', 'ricc', if_lyapunov)
   if (main_run) then
      write(home,'(A,A,I3.3,A,A,A)') trim(home_base), 'Tend', int(Tend), '_', eq, '/'
   else
      home = home_base
   end if

   ! Setup logging
   call logger_setup(logfile=trim(home)//'lightkrylov.log', log_level=error_level, log_stdout=.false., log_timestamp=.true.)

   ! Initialize timers for LightKrylov and LightROM
   call initialize_timers()
   call global_lightROM_timer%add_timer('DLRA Ginzburg-Landau example', start=.true.)
   
   ! Enumerate timers to check proper initialization
   call enumerate_timers()

   call print_test_info(if_lyapunov, if_adj, main_run, run_fixed_rank_test, run_rank_adaptive_test, run_eigenvalue_test)
   call print_test_header(if_lyapunov, if_adj, eq, refid, rk_X0, Xref)
   
   ! Initialize propagator
   print '(4X,A)', 'Initialize exponential propagator'
   prop = exponential_prop(1.0_dp)
   
   ! Initialize LTI system
   A = GL_operator()
   ! Initialize mesh and system parameters A, B, CT
   print '(4X,A)', 'Initialize LTI system (A, prop, B, CT, _)'
   LTI = lti_system()
   call LTI%initialize_lti_system(A, prop, B, CT)
      
   ! Define initial condition
   allocate(U0(rk_X0), source=B(1)); call zero_basis(U0)
   allocate(S0(rk_X0,rk_X0)); S0 = 0.0_dp
   print '(4X,A)', 'Define initial condition'
   call generate_random_initial_condition(U0, S0, rk_X0)
   call reconstruct_solution(X_out, U0, S0)
   print *, ''
   print '(A,F16.12)', '  |  X_0  |/N = ', norm2(X_out)/N
   print '(A,F16.12)', '  | res_0 |/N = ', norm2(CALE(X_out, if_adj))/N
   print *, ''
   call print_svdvals(X_out, ' X0 ', 60)
   
   if (main_run) then
      Tend = Tend_long
      rkv  = [ 10, 14, 18 ]
      dtv  = logspace(-3.0_dp, 0.0_dp, 4, 10)
      dtv  = dtv(size(dtv):1:-1) ! reverse vector
      tolv = [ 1e-2_dp, 1e-6_dp, 1e-10_dp ]
      if (if_lyapunov) then
         TOv = [ 1, 2 ]
      else
         TOv = [ 1 ]
      end if
   else
      Tend = Tend_short
      dtv  = logspace(-4.0_dp, -2.0_dp, 3, 10)
      dtv  = dtv(size(dtv):1:-1) ! reverse vector
      tolv = [ 1e-2_dp ]
      if (if_lyapunov) then
         rkv = [ 8, 12, 16 ]
         TOv = [ 1, 2 ]
      else
         rkv = [ 6, 10, 14 ]
         TOv = [ 1 ]
      end if
   end if
   nprint = 60
  
   !
   !        RK long time horizon
   !
   call integrate_RK(eq, LTI, Xref, Xref_RK, U0, S0, Tend, if_adj, home_base, main_run)
   
   !
   !        rank-adaptive DLRA long time horizon
   !
   if (run_fixed_rank_test) then
      call integrate_DLRA_fixed(eq, LTI, Xref, Xref_RK, U0, S0, Tend, dtv, rkv, TOv, nprint, if_adj, home, main_run)
   end if
   
   !
   !        rank-adaptive DLRA long time horizon
   !
   if (run_rank_adaptive_test) then
      call integrate_DLRA_adaptive(eq, LTI, Xref, Xref_RK, U0, S0, Tend, dtv, TOv, tolv, nprint, if_adj, home, main_run)
   end if

   if (if_lyapunov) then
      tmr_name = 'BPOD'
      call global_lightROM_timer%add_timer(tmr_name, start=.true.)
      Tendv = [ 10.0_dp, 50.0_dp, 100.0_dp ]
      dtv  = logspace(-1.0_dp, 0.0_dp, 3, 10)
      dtv  = dtv(size(dtv):1:-1) ! reverse vector
      print '(1X,A)', 'POD:'
      do it = 1, size(tolv)
         T_POD = tolv(it)
         print '(3X,A,F12.6)', '   Tend:', T_POD
         do irep = 1, size(dtv)
            tau = dtv(irep)
            prop = exponential_prop(tau)
            if (if_adj) then
               allocate(X0(rk_C), source=CT)
            else
               allocate(X0(rk_B), source=B)
            end if
            call Proper_Orthogonal_Decomposition(svals, prop, X0, tau, T_POD, if_adj)
            nprint = min(8, size(svals))
            call print_svdvals(svals, 'XTX', nprint)
            deallocate(X0)
         end do
      end do
      call global_lightROM_timer%stop(tmr_name)
   end if

   if (run_eigenvalue_test .and. .not. if_lyapunov) then
      dtv  = logspace(-4.0_dp, 0.0_dp, 5, 10)
      dtv  = dtv(size(dtv):1:-1) ! reverse vector
      tolv = [ 1e-2_dp, 1e-6_dp, 1e-10_dp ]
      torder = 1
      
      ! eig A
      tmr_name =  'eig A'
      fname    = trim(home)//"spectrum_A.npy"
      call eigenvalue_analysis(prop, tmr_name, fname)
      print *, ''
      
      ! eig A - BK
      X = LR_state()
      if (if_adj) then
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
      call prop_control%init(X, B, Rinv, adjoint=if_adj, enable_control=.true.)
      
      call eigenvalue_analysis(prop_control, tmr_name, fname)

      ! deallocate and clean
      deallocate(rkintegrator); deallocate(prop_control)
      
      print *, ''
      do j = 1, size(tolv)
         tol = tolv(j)
         do k = 1, size(dtv)
            tau = dtv(k)
            note = merge('Padj', 'Pdir', if_adj)
            fbase = make_filename(home, 'DLRA_ADAPT', eq, trim(note), rk, torder, tau, Tend, tol)
            exist_file = exist_X_file(fbase)
            if (exist_file) then
               ! load X state
               call load_X_from_file(X, meta, fbase, U0)
               ! recreate integrators
               rkintegrator = rks54_class_with_control()
               prop_control = exponential_prop_with_control(1.0_dp, prop=rkintegrator)
               ! initialize integrators
               call prop_control%init(X, B, Rinv, adjoint=if_adj, enable_control=.true.)
               
               call make_labels(Tstr, taustr, tolstr, Tend, tau, tol)
               note = merge('A-LC', 'A-BK', if_adj)
               fname = trim(home)//'spectrum_'//trim(note)//'_Tend'//trim(Tstr)//'_tau'//trim(taustr)//'_tol'//trim(tolstr)//'.npy'
               write(tmr_name,'(*(A))')  'eig ', trim(note), ':   Tend= ', trim(Tstr), '   tau= ', trim(taustr), '   tol= ', trim(tolstr)
               
               ! eigendecomposition
               call eigenvalue_analysis(prop_control, tmr_name, fname)
               
               ! deallocate and clean
               deallocate(rkintegrator); deallocate(prop_control)
            end if
         end do
         write(*,*)
      end do
   end if

   ! Compute and print timer summary
   call finalize_timers()

   return
end program demo