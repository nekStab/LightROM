program demo
   ! Standard Library.
   use stdlib_optval, only : optval 
   use stdlib_linalg, only : eye, diag
   use stdlib_math, only : all_close, logspace
   use stdlib_io_npy, only : save_npy, load_npy
   use stdlib_logger, only : information_level, warning_level, debug_level, error_level, none_level
    ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov, only : wp => dp
   use LightKrylov_Constants
   use LightKrylov_Logger
   use LightKrylov_AbstractVectors
   use LightKrylov_ExpmLib
   use LightKrylov_Utils
   ! LightROM
   use LightROM_AbstractLTIsystems
   use LightROM_Utils
   use LightROM_LyapunovSolvers
   use LightROM_LyapunovUtils
   ! GInzburg-Landau
   use Ginzburg_Landau_Base
   use Ginzburg_Landau_Operators
   use Ginzburg_Landau_Utils
   use Ginzburg_Landau_Tests
   implicit none

   ! DLRA
   logical, parameter :: verb  = .true.
   !
   logical :: run_test
   !
   character(len=128)      :: oname
   character(len=128)      :: onameU
   character(len=128)      :: onameS
   character(len=128)      :: testinfo
   ! rk_B & rk_C are set in ginzburg_landau_base.f90

   integer  :: nrk, ntau, rk,  torder
   real(wp) :: tau, Tend, Ttot
   ! vector of dt values
   real(wp), allocatable :: tauv(:), tolv(:), taucv(:)
   ! vector of rank values
   integer, allocatable :: rkv(:), TOv(:)

   ! Exponential propagator (RKlib).
   type(GL_operator),      allocatable       :: A
   type(exponential_prop), allocatable       :: prop

   ! LTI system
   type(lti_system)                          :: LTI

   ! Initial condition
   type(state_vector)                        :: U0(rk_X0)
   real(wp)                                  :: S0(rk_X0,rk_X0)
   ! matrix
   real(wp)                                  :: U0_in(2*nx, rkmax)
   
   ! OUTPUT
   real(wp)                                  :: U_out(2*nx,rkmax)
   real(wp)                                  :: X_out(2*nx,2*nx)
   real(wp)                                  :: lagsvd(rkmax)

   ! Information flag.
   integer                                   :: info

   ! Counters
   integer                                   :: i, j, k, irep, nrep, istep, nsteps
   integer,                allocatable       :: perm(:)
   real(wp)                                  :: Tmax, etime

   logical        :: ifsave, ifload, ifverb, iflogs, ifrk, ifconv

   call logger_setup()
   
   !----------------------------------
   !-----     INITIALIZATION     -----
   !----------------------------------

   ! Initialize mesh and system parameters A, B, CT
   if (verb .and. io_rank()) write(*,*) 'Initialize parameters'
   call initialize_parameters()

   ! Initialize propagator
   if (verb .and. io_rank()) write(*,*) 'Initialize exponential propagator'
   prop = exponential_prop(1.0_wp)

   ! Initialize LTI system
   A = GL_operator()
   if (verb .and. io_rank()) write(*,*) 'Initialize LTI system (A, prop, B, CT, _)'
   LTI = lti_system()
   call LTI%initialize_lti_system(A, prop, B, CT)

   ! Define initial condition of the form X0 + U0 @ S0 @ U0.T SPD 
   if (verb .and. io_rank()) write(*,*) '    Define initial condition'
   call generate_random_initial_condition(U0, S0, rk_X0)
   call get_state(U_out(:,:rk_X0), U0)
   
   if (verb .and. io_rank()) write(*,*) '    Start tests'
   if (verb .and. io_rank()) write(*,*) '    Output folder:', basepath

   !----------------------------------
   !
   ! DLRA CONVERGENCE TEST FOR LYAPUNOV EQUATION
   !
   !----------------------------------

   run_test = .false.
   if (run_test) then
      nrk  = 8; allocate(rkv(1:nrk));   rkv  = (/ 2, 6, 10, 14, 20, 40, 80, 128, 256 /)
      ntau = 3; allocate(tauv(1:ntau)); tauv = logspace(-4.0, -3.0, ntau)
      allocate(TOv(2)); TOv = (/ 1, 2 /)
      Tend = 0.01_wp
      ! run DLRA
      ifsave = .true. ! save X_rk to disk (LightROM/local)
      ifverb = .true. ! verbosity
      call run_lyap_convergence_test(LTI, U0, S0, Tend, tauv, rkv, TOv, ifverb)
      deallocate(rkv)
      deallocate(tauv)
      deallocate(TOv)
   end if 

   !----------------------------------
   !
   ! DLRA TEST FOR LYAPUNOV EQUATION
   !
   !----------------------------------

   run_test = .false.
   if (run_test) then
      !nrk  = 6; allocate(rkv(1:nrk));   rkv  = (/ 6, 10, 12, 14, 20, 40 /)
      !ntau = 5; allocate(tauv(1:ntau)); tauv = (/ 1.0, 0.1, 0.01, 0.001, 0.0001 /)
      !allocate(TOv(2)); TOv = (/ 1, 2 /)
      nrk  = 1; allocate(rkv(1:nrk));   rkv  = (/ 12 /)
      ntau = 1; allocate(tauv(1:ntau)); tauv = (/ 0.1 /)
      allocate(TOv(1)); TOv = (/ 1 /)
      Tend = 1.0_wp
      nrep = 60
      ! run DLRA
      ifsave = .false. ! save X and Y matrices to disk (LightROM/local)
      ifverb = .false. ! verbosity
      iflogs = .false. ! write logs with convergence and signular value evolution
      call run_DLRA_lyapunov_test(LTI, U0, S0, rkv, tauv, TOv, Tend, nrep, ifsave, ifverb, iflogs)
      deallocate(rkv)
      deallocate(tauv)
      deallocate(TOv)
   end if 

   !----------------------------------
   !
   ! DLRA TEST FOR BALANCING TRANSFORMATION
   !
   !----------------------------------

   run_test = .false.
   if (run_test) then
      ! Set parameters
      rk     = 14
      tau    = 0.1_wp
      torder = 2
      ! integration time
      Tmax   = 60.0_wp
      nrep   = 12
      Tend   = Tmax/nrep
      nsteps = nint(Tend/tau)
      ! run DLRA
      ifsave = .true.  ! save X and Y matrices to disk (LightROM/local)
      ifload = .false. ! read X and Y matrices from disk (LightROM/local)
      ifverb = .true.  ! verbosity
      iflogs = .true.  ! write logs with convergence and signular value evolution
      call run_BT_test(LTI, U0, S0, rk, tau, torder, Tmax, nrep, ifsave, ifload, ifverb, iflogs)
   end if

   !----------------------------------
   !
   ! TEST COMPARING THE SPEED OF RK vs KRYLOV EXPONTNTIAL INTEGRATORS
   !
   !----------------------------------

   run_test = .false.
   if (run_test) then
      ntau = 20; allocate(tauv(1:ntau)); tauv = logspace(-5.0_wp,0.0_wp,ntau)
      torder = 1
      call run_kexpm_test(LTI%A, LTI%prop, U0(1), tauv, torder, 1000)
   end if

   !----------------------------------
   !
   ! TEST COMPARING THE SPEED OF RK vs KRYLOV EXPONTNTIAL INTEGRATORS
   !
   !----------------------------------

   run_test = .false.
   if (run_test) then
      ntau = 5; allocate(tauv(1:ntau)); tauv = logspace(-5.0_wp,0.0_wp,ntau)
      torder = 1
      call run_kexpm_var_dt_test(LTI%A, LTI%prop, U0(1), tauv, torder, 100)
   end if

   !----------------------------------
   !
   ! DLRA TEST FOR RICCATI EQUATION
   !
   !----------------------------------

   run_test = .false.
   if (run_test) then
      nrk  = 6; allocate(rkv(1:nrk));   rkv  = (/ 6, 10, 12, 14, 20, 40 /)
      ntau = 5; allocate(tauv(1:ntau)); tauv = (/ 1.0, 0.1, 0.01, 0.001, 0.0001 /)
      Tend = 1.0_wp
      nrep = 60
      ! run DLRA
      ifsave = .true. ! save X and Y matrices to disk (LightROM/local)
      ifverb = .true. ! verbosity
      iflogs = .true. ! write logs with convergence and signular value evolution
      call run_DLRA_riccati_test(LTI, U0, S0, Qc, Rinv, &
                                 & rkv, tauv, Tend, nrep, & 
                                 & ifsave, ifverb, iflogs)
      deallocate(rkv)
      deallocate(tauv)
   end if 

   !----------------------------------
   !
   ! RANK-ADAPTIVE DLRA LYAPUNOV
   !
   !----------------------------------

   run_test = .true.
   if (run_test) then
      rkv  = (/ 6 /)
      tauv = (/ 0.1 /)
      tolv = (/ 1e-6 /)
      TOv = (/ 1 /)
      Tend = 60.0_wp
      ! run DLRA
      ifsave = .false. ! save X and Y matrices to disk (LightROM/local)
      ifverb = .false. ! verbosity
      ifrk   = .true.  ! use Runge-Kutto (or kexpm)
      ifconv = .false. ! convergence run? 
      call logger%configure(level=error_level)
      call run_DLRA_rank_adaptive_test(LTI, U0, S0, rkv, tauv, TOv, tolv, Tend, 'dlra_test', ifsave, ifverb, ifrk, ifconv)
   end if 

   !----------------------------------
   !
   ! RANK-ADAPTIVE DLRA LYAPUNOV with convergence run
   !
   !----------------------------------

   run_test = .false.
   if (run_test) then
      call logger%configure(level=warning_level)
      ifsave = .true. ! save X and Y matrices to disk (LightROM/local)
      ifverb = .true. ! verbosity
      iflogs = .true. ! write logs with convergence and signular value evolution
      Tend = 50.0_wp
      ! 01
      testinfo = '01_rk_tauv'
      rkv  = (/ 4, 6, 8, 10 /)
      tauv = (/ 0.001, 0.01, 0.1, 1.0 /)
      tolv = (/ 1e-6 /)
      TOv  = (/ 1 /)
      ! run DLRA
      call run_DLRA_rank_adaptive_test(LTI, U0, S0, rkv, tauv, TOv, tolv, Tend, testinfo, ifsave, ifverb, .false.)
      ! exclude small dt that do not converge to right rank
      !rkv  = (/ 10 /) 
      rkv = (/ 4, 6, 8, 10 /)
      !tauv = (/ 1.0 /) 
      tauv = (/ 0.1, 1.0 /)
      call run_DLRA_rank_adaptive_test(LTI, U0, S0, rkv, tauv, TOv, tolv, Tend, testinfo, ifsave, ifverb, .true.)
      ! 02
      testinfo = '02_rk_tauv_tol'
      rkv  = (/ 6, 10 /)
      tauv = (/ 0.001, 0.01, 0.1, 1.0 /)
      tolv = (/ 1e-8, 1e-12 /)
      TOv  = (/ 1 /)
      ! run DLRA
      call run_DLRA_rank_adaptive_test(LTI, U0, S0, rkv, tauv, TOv, tolv, Tend, testinfo, ifsave, ifverb, .false.)
      ! exclude small dt that do not converge to right rank
      !rkv  = (/ 10 /) 
      rkv = (/ 6, 10 /)
      !tauv = (/ 1.0 /) 
      tauv = (/ 0.1, 1.0 /)
      tolv = (/ 1e-12 /)
      call run_DLRA_rank_adaptive_test(LTI, U0, S0, rkv, tauv, TOv, tolv, Tend, testinfo, ifsave, ifverb, .true.)

      ! 03
      testinfo = '03_TO2'
      rkv  = (/ 10 /)
      tauv = (/ 0.001, 0.01, 0.1, 1.0 /)
      tolv = (/ 1e-6 /)
      TOv  = (/ 2 /)
      ! run DLRA
      call run_DLRA_rank_adaptive_test(LTI, U0, S0, rkv, tauv, TOv, tolv, Tend, testinfo, ifsave, ifverb, .false.)
      ! exclude small dt that do not converge to right rank
      tauv = (/ 0.1, 1.0 /)
      nrep = 50
      call run_DLRA_rank_adaptive_test(LTI, U0, S0, rkv, tauv, TOv, tolv, Tend, testinfo, ifsave, ifverb, .true.)
   end if


   return
end program demo
