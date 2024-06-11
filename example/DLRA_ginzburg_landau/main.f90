program demo
   ! Standard Library.
   use stdlib_optval, only : optval 
   use stdlib_linalg, only : eye, diag
   use stdlib_math, only : all_close, logspace
   use stdlib_io_npy, only : save_npy, load_npy
    ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov, only : wp => dp
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
   character*128      :: oname
   character*128      :: onameU
   character*128      :: onameS
   ! rk_B & rk_C are set in ginzburg_landau_base.f90

   integer  :: nrk, ntau, rk,  torder
   real(wp) :: tau, Tend, Ttot
   ! vector of dt values
   real(wp), allocatable :: tauv(:)
   ! vector of rank values
   integer, allocatable :: rkv(:), TOv(:)

   ! Exponential propagator (RKlib).
   type(GL_operator),      allocatable       :: A
   type(exponential_prop), allocatable       :: prop

   ! LTI system
   type(lti_system)                          :: LTI

   ! Initial condition
   type(state_vector)                        :: U0(1:rkmax)
   real(wp)                                  :: S0(rkmax,rkmax)
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

   logical        :: ifsave, ifload, ifverb, iflogs
   
   call logger%configure(level=none_level)
   
   !----------------------------------
   !-----     INITIALIZATION     -----
   !----------------------------------

   ! Initialize mesh and system parameters A, B, CT
   if (verb) write(*,*) 'Initialize parameters'
   call initialize_parameters()

   ! Initialize propagator
   if (verb) write(*,*) 'Initialize exponential propagator'
   prop = exponential_prop(1.0_wp)

   ! Initialize LTI system
   A = GL_operator()
   if (verb) write(*,*) 'Initialize LTI system (A, prop, B, CT, _)'
   LTI = lti_system()
   call LTI%initialize_lti_system(A, prop, B, CT)

   ! Define initial condition of the form X0 + U0 @ S0 @ U0.T SPD 
   if (verb) write(*,*) '    Define initial condition'
   call generate_random_initial_condition(U0, S0, rk_X0)
   call get_state(U_out, U0)

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

   run_test = .true.
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
      ifsave = .true. ! save X and Y matrices to disk (LightROM/local)
      ifverb = .true. ! verbosity
      iflogs = .true. ! write logs with convergence and signular value evolution
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
!
   return
end program demo