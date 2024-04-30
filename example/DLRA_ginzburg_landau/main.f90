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
   logical, parameter :: if_run_DLRA_test = .true.
   logical, parameter :: if_run_BT_test = .false.
   logical, parameter :: if_run_conv_metrics_test = .false.
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
   integer, allocatable :: rkv(:)

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

   if (if_run_DLRA_test) then
      nrk  = 6; allocate(rkv(1:nrk));   rkv  = (/ 2, 6, 10, 14, 20, 40 /)
      ntau = 5; allocate(tauv(1:ntau)); tauv = (/ 0.5, 0.1, 0.01, 0.001, 0.0001 /)
      Tend = 1.0_wp
      nrep = 60
      ! run DLRA
      ifsave = .true. ! save X and Y matrices to disk (LightROM/local)
      ifverb = .true. ! verbosity
      iflogs = .true. ! write logs with convergence and signular value evolution
      call run_DLRA_test(LTI, U0, S0, rkv, tauv, Tend, nrep, ifsave, ifverb, iflogs)
   end if 

   if (if_run_BT_test) then
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
!
   return
end program demo