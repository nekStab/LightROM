program demo
   use LightKrylov
   use LightKrylov_expmlib
   use LightKrylov_utils

   use LightROM_AbstractLTIsystems
   use LightROM_utils

   use LightROM_LyapunovSolvers
   use LightROM_LyapunovUtils

   use Ginzburg_Landau_Base
   use Ginzburg_Landau_Operators
   use Ginzburg_Landau_Utils

   use stdlib_optval, only : optval 
   use stdlib_linalg, only : eye, diag
   use stdlib_math, only : all_close, logspace
   use stdlib_io_npy, only : save_npy, load_npy

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
   real(kind=wp) :: tau, Tend, Ttot
   ! vector of dt values
   real(kind=wp), allocatable :: tauv(:)
   ! vector of rank values
   integer, allocatable :: rkv(:), TOv(:)

   ! Exponential propagator (RKlib).
   type(GL_operator),      allocatable       :: A
   type(exponential_prop), allocatable       :: prop

   ! LTI system
   type(lti_system)                          :: LTI

   ! Initial condition
   real(kind=wp)                             :: U0_in(2*nx, rkmax)
   real(kind=wp)                             :: S0(rkmax,rkmax)
   type(state_vector), allocatable           :: U0(:)
   type(state_vector), allocatable           :: Utmp(:)
   
   ! OUTPUT
   real(kind=wp)                             :: U_out(2*nx,rkmax)
   real(kind=wp)                             :: X_out(2*nx,2*nx)
   real(kind=wp)                             :: lagsvd(rkmax)

   ! Riccati
   real(kind=wp)                             :: Qc(rk_c,rk_c)
   real(kind=wp)                             :: Rinv(rk_b, rk_b)
   real(kind=wp)                             :: sqrtw(2*nx)

   ! Information flag.
   integer                                   :: info

   ! Counters
   integer                                   :: i, j, k, irep, nrep, istep, nsteps
   real(kind=wp)                             :: Tmax, etime
!
   ! SVD
   real(kind=wp)  :: U_svd(2*nx,2*nx)
   real(kind=wp)  :: S_svd(rkmax)
   real(kind=wp)  :: V_svd(rkmax,rkmax)

   logical        :: ifsave, ifload, ifverb, iflogs
   
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

   ! Define initial condition
   if (verb) write(*,*) 'Define initial condition'
   allocate(U0(1:rk_X0))
   call init_rand(U0, .false.)
   call get_state(U0_in(:,1:rk_X0), U0)
   ! Compute SVD to get low-rank representation
   call svd(U0_in(:,1:rk_X0), U_svd(:,1:2*nx), S_svd(1:rk_X0), V_svd(1:rk_X0,1:rk_X0))
   S0 = 0.0_wp
   S0(1:rk_X0,1:rk_X0) = diag(S_svd(1:rk_X0))
   call set_state(U0, U_svd(:,1:rk_X0))

   !call get_state(U_out(:,1:1), LTI%B)
   !sqrtw = sqrt(weight)
   !U_out(:,1) = sqrtw*U_out(:,1)
   !X_out = 0.0_wp
   !X_out = matmul(U_out(:,1:1), transpose(U_out(:,1:1)))
   !call print_mat(2*nx, 2*nx, X_out)
   !print *, maxval(X_out)
   !print *, weight(1:10)
   !STOP 1

   !----------------------------------
   !
   ! DLRA TEST FOR LYAPUNOV EQUATION
   !
   !----------------------------------

   run_test = .true.
   if (run_test) then
      nrk  = 6; allocate(rkv(1:nrk));   rkv  = (/ 2, 6, 10, 14, 20, 40 /)
      ntau = 4; allocate(tauv(1:ntau)); tauv = (/ 1.0, 0.1, 0.01, 0.001 /)
      !ntau = 1; allocate(tauv(1:ntau)); tauv = (/ 1.0 /)
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

   STOP 1
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
      Tmax   = 50.0_wp
      nrep   = 10
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
      nrk  = 6; allocate(rkv(1:nrk));   rkv  = (/ 2, 6, 10, 14, 20, 40 /)
      ntau = 2; allocate(tauv(1:ntau)); tauv = (/ 1.0, 0.1 /)
      Tend = 1.0_wp
      nrep = 60
      Qc   = eye(rk_c)
      Rinv = eye(rk_b)
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