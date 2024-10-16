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
   use Ginzburg_Landau_RK
   implicit none

   ! DLRA
   logical, parameter :: verb  = .true.
   !
   logical :: run_test
   !
   character*128      :: onameU, onameS, oname
   ! rk_B & rk_C are set in ginzburg_landau_base.f90

   integer  :: nrk, ntau, rk, torder, nsteps, iref
   real(wp) :: tau, Tend, T_RK
   ! vector of dt values
   real(wp), allocatable :: tauv(:), tolv(:)
   ! vector of rank values
   integer, allocatable :: rkv(:), TOv(:)

   ! Exponential propagator (RKlib).
   type(GL_operator),      allocatable       :: A
   type(exponential_prop), allocatable       :: prop

   ! LTI system
   type(lti_system)                          :: LTI

   ! Initial condition
   type(state_vector),     allocatable       :: U0(:)
   real(wp),               allocatable       :: S0(:,:)
   ! matrix
   real(wp)                                  :: U0_in(2*nx, rkmax)
   
   ! OUTPUT
   real(wp)                                  :: U_out(2*nx,rkmax)
   real(wp)                                  :: X_out(2*nx,2*nx)
   real(wp)                                  :: lagsvd(rkmax)
   real(wp)                                  :: res_flat(N**2)

   ! Reference solutions (BS & RK)
   real(wp)                                  :: Xref_BS(N,N)
   real(wp)                                  :: Xref_RK(N,N)

   ! IO
   real(wp),           allocatable           :: U_load(:,:)

   ! Information flag.
   integer                                   :: info

   ! Counters
   integer                                   :: i, j, k, irep, nrep, istep
   real(wp)                                  :: etime

   real(wp)                           :: U_svd(2*nx,2*nx)
   real(wp)                           :: S_svd(2*nx)
   real(wp)                           :: V_svd(2*nx,2*nx)

   logical        :: ifsave, ifload, iflogs

   ! Exact comparison
   type(dlra_opts)                              :: opts
   type(LR_state),                allocatable   :: X
   class(abstract_vector_rdp),             allocatable   :: exptAU(:)    ! scratch basis
   real(wp),                               allocatable   :: R(:,:)    ! QR coefficient matrix
   integer,                                allocatable   :: perm(:)   ! Permutation vector
   real(wp),                               allocatable   :: wrk(:,:)   

   !call logger%configure(level=information_level, time_stamp=.false.)
   call logger%configure(level=error_level, time_stamp=.false.)
   
   !----------------------------------
   !-----     INITIALIZATION     -----
   !----------------------------------

   ! Initialize mesh and system parameters A, B, CT
   if (verb) print *, 'Initialize parameters'
   call initialize_parameters()

   ! Initialize propagator
   if (verb) print *, 'Initialize exponential propagator'
   prop = exponential_prop(1.0_wp)

   ! Initialize LTI system
   A = GL_operator()
   if (verb) print *, 'Initialize LTI system (A, prop, B, CT, _)'
   LTI = lti_system()
   call LTI%initialize_lti_system(A, prop, B, CT)

   print *, ''
   print *, 'Check residual computation with Bartels-Stuart solution (python):'
   oname = trim(basepath)//"BS/CGL_Lyapunov_Controllability_Xref_BS_W.npy"
   call load_data(oname, U_load)
   Xref_BS = U_load
   
   call CALE(res_flat, reshape(Xref_BS, shape(res_flat)), BBTW_flat, .false.)
   print *, ''
   print '(A,F18.12)', '   ||  X_BS  ||_2/N = ', norm2(Xref_BS)/N
   print '(A,F18.12)', '   || res_BS ||_2/N = ', norm2(res_flat)/N
   print *, ''
   
   ! Define initial condition of the form X0 = U0 @ S0 @ U0.T (SPD) or load from file
   ifload = .true.
   onameU = 'local/CGL_Nx256_U0_rk_X0_10.npy'
   onameS = 'local/CGL_Nx256_S0_rk_X0_10.npy'
   allocate(U0(rk_X0), source=B(1)); call zero_basis(U0)
   allocate(S0(rk_X0,rk_X0)); S0 = 0.0_wp
   if (ifload) then
      if (verb) print *, 'Load initial condition'
      call load_data(onameU, U_load)
      call set_state(U0, U_load(:,:rk_X0)/sqrt(weight(1)))
      call load_data(onameS, U_load)
      S0(:rk_X0,:rk_X0) = U_load(:rk_X0,:rk_X0)
      call get_state(U_out(:,:rk_X0), U0)
   else
      if (verb) print *, 'Define initial condition'
      call generate_random_initial_condition(U0, S0, rk_X0)
      call get_state(U_out(:,:rk_X0), U0)
      call save_data(onameU, U_out(:,:rk_X0))
      call save_data(onameS, S0(:rk_X0,:rk_X0))
   end if
   X_out = matmul(matmul( U_out(:,:rk_X0), matmul( S0(:rk_X0,:rk_X0), transpose(U_out(:,:rk_X0)) ) ), weight_mat)
   call CALE(res_flat, reshape(X_out, shape(res_flat)), BBTW_flat, .false.)
   print *, ''
   print *, 'Cross-check initial condition:'
   print *, ''
   print '(A,F18.12)', '    ||  X_0  ||_2/N = ', norm2(X_out)/N
   print '(A,F18.12)', '    || res_0 ||_2/N = ', norm2(res_flat)/N
   ! compute svd
   call svd(X_out, S_svd, U_svd, V_svd)
   print *, 'SVD:'
   print '(3X,*(F18.12,1X))', ( S_svd(i), i = 1, 5)
   print '(3X,*(F18.12,1X))', ( S_svd(i), i = 6, 10)
   
   ! Run RK integrator for the Lyapunov equation
   T_RK   = 1.0_wp
   nsteps = 10
   iref   = 1
   ifload = .true.
   call run_lyap_reference_RK(LTI, Xref_BS, Xref_RK, U0, S0, T_RK, nsteps, iref, ifload)

   call CALE(res_flat, reshape(Xref_RK, shape(res_flat)), BBTW_flat, .false.)
   print *, ''
   print '(A,F18.12)', '   ||  X_RK  ||_2/N = ', norm2(Xref_BS)/N
   print '(A,F18.12)', '   || res_RK ||_2/N = ', norm2(res_flat)/N
   print *, ''
   
   ! basic settings
   Tend = T_RK/nsteps*iref
   ifsave = .false. ! save data to disk (LightROM/local)

   ! Testing
   rk = 4
   tau = 0.01_wp
   opts = dlra_opts(chktime=1.0_wp, inc_tol=atol_dp, if_rank_adaptive=.false.)
   X = LR_state()

   ! Initialize low-rank representation with rank rk
   call X%initialize_LR_state(U0, S0, rk, rkmax, .false.)
   nsteps = nint(Tend/tau)

   ! first order
   ! M Step
   rk = X%rk
   allocate(R(rk,rk));   R    = 0.0_wp 
   allocate(perm(rk));   perm = 0
   allocate(wrk(rk,rk)); wrk  = 0.0_wp
   allocate(exptAU(rk), source=X%U(1))
   ! Apply propagator to initial basis using k_exptA
   call zero_basis(exptAU)
   do i = 1, rk
      call k_exptA_rdp(exptAU(i), A, X%U(i), tau, info, .false.)
   end do
   call get_state(U_out(:,:rk), exptAU(:rk))
   call save_data(trim(basepath)//'TEST/exptA_kexpm.npy', U_out(:,:rk))
   ! Apply propagator to initial basis using RKlib
   call zero_basis(exptAU)
   do i = 1, rk
      call exptA(exptAU(i), A, X%U(i), tau, info, .false.)
      call X%U(i)%axpby(0.0_wp, exptAU, 1.0_wp) ! overwrite old solution
   end do
   call get_state(U_out(:,:rk), exptAU(:rk))
   call save_data(trim(basepath)//'TEST/exptA_RK.npy', U_out(:,:rk))
   ! Reorthonormalize in-place without pivoting
   call qr(X%U(:rk), R, info)
   call get_state(U_out(:,:rk), X%U(:rk))
   call save_data(trim(basepath)//'TEST/qr_Q.npy', U_out(:,:rk))
   call save_data(trim(basepath)//'TEST/qr_R.npy', R)
   ! reset data
   do i = 1, rk
      call X%U(i)%axpby(0.0_wp, exptAU, 1.0_wp) ! overwrite old solution
   end do
   ! Reorthonormalize in-place with pivoting
   call qr(X%U(:rk), R, perm, info)
   call apply_inverse_permutation_matrix(R, perm)
   call get_state(U_out(:,:rk), X%U(:rk))
   call save_data(trim(basepath)//'TEST/qrp_Q.npy', U_out(:,:rk))
   call save_data(trim(basepath)//'TEST/qrp_R.npy', R)

   STOP 6
   ! DLRA with fixed rank
   rkv = (/ 4, 20 /)
   tauv = logspace(-4.0, -1.0, 4, 10)
   !tauv = tauv(size(tauv):1:-1) ! reverse vector
   TOv  = (/ 1, 2 /)
   call run_lyap_DLRA_test(LTI, Xref_BS, Xref_RK, U0, S0, Tend, tauv, rkv, TOv, ifsave)
   
   ! DLRA with adaptive rank
   !rkv = (/ 4, 20 /) !(/ 20 /) !(/ 4, 8, 10, 20, 40/)
   !tauv = (/ 0.001 /) !tauv = logspace(-3.0_wp, 0.0_wp, 4, 10) !(/ 0.1 /) !logspace(-3.0, 0.0, 7, 10)
   !!tauv = tauv(size(tauv):1:-1) ! reverse vector
   !TOv  = (/ 1 /)  !(/ 1, 2 /)
   !tolv = (/ 1e-2_wp, 1e-6_wp, 1e-10_wp /)
   !call run_lyap_DLRArk_test(LTI, Xref_BS, Xref_RK, U0, S0, Tend, tauv, rkv, TOv, tolv, ifsave)
   !STOP 6
   
   !----------------------------------
   !
   ! DLRA CONVERGENCE TEST FOR LYAPUNOV EQUATION
   !
   !----------------------------------

   !run_test = .false.
   !if (run_test) then
   !   nrk  = 8; allocate(rkv(1:nrk));   rkv  = (/ 2, 6, 10, 14, 20, 40, 80, 128, 256 /)
   !   ntau = 3; allocate(tauv(1:ntau)); tauv = logspace(-4.0, -3.0, ntau)
   !   allocate(TOv(2)); TOv = (/ 1, 2 /)
   !   Tend = 0.01_wp
   !   ! run DLRA
   !   ifsave = .true. ! save X_rk to disk (LightROM/local)
   !   ifverb = .true. ! verbosity
   !   call run_lyap_convergence_test(LTI, U0, S0, Tend, tauv, rkv, TOv, ifverb)
   !   deallocate(rkv)
   !   deallocate(tauv)
   !   deallocate(TOv)
   !end if 

   !----------------------------------
   !
   ! DLRA TEST FOR LYAPUNOV EQUATION
   !
   !----------------------------------

   !run_test = .false.
   !if (run_test) then
   !   !nrk  = 6; allocate(rkv(1:nrk));   rkv  = (/ 6, 10, 12, 14, 20, 40 /)
   !   !ntau = 5; allocate(tauv(1:ntau)); tauv = (/ 1.0, 0.1, 0.01, 0.001, 0.0001 /)
   !   !allocate(TOv(2)); TOv = (/ 1, 2 /)
   !   nrk  = 1; allocate(rkv(1:nrk));   rkv  = (/ 12 /)
   !   ntau = 1; allocate(tauv(1:ntau)); tauv = (/ 0.1 /)
   !   allocate(TOv(1)); TOv = (/ 1 /)
   !   Tend = 1.0_wp
   !   nrep = 60
   !   ! run DLRA
   !   ifsave = .false. ! save X and Y matrices to disk (LightROM/local)
   !   ifverb = .false. ! verbosity
   !   iflogs = .false. ! write logs with convergence and signular value evolution
   !   call run_DLRA_lyapunov_test(LTI, U0, S0, rkv, tauv, TOv, Tend, nrep, ifsave, ifverb, iflogs)
   !   deallocate(rkv)
   !   deallocate(tauv)
   !   deallocate(TOv)
   !end if 

   !----------------------------------
   !
   ! DLRA TEST FOR BALANCING TRANSFORMATION
   !
   !----------------------------------

   !run_test = .false.
   !if (run_test) then
   !   ! Set parameters
   !   rk     = 14
   !   tau    = 0.1_wp
   !   torder = 2
   !   ! integration time
   !   Tmax   = 60.0_wp
   !   nrep   = 12
   !   Tend   = Tmax/nrep
   !   nsteps = nint(Tend/tau)
   !   ! run DLRA
   !   ifsave = .true.  ! save X and Y matrices to disk (LightROM/local)
   !   ifload = .false. ! read X and Y matrices from disk (LightROM/local)
   !   ifverb = .true.  ! verbosity
   !   iflogs = .true.  ! write logs with convergence and signular value evolution
   !   call run_BT_test(LTI, U0, S0, rk, tau, torder, Tmax, nrep, ifsave, ifload, ifverb, iflogs)
   !end if

   !----------------------------------
   !
   ! TEST COMPARING THE SPEED OF RK vs KRYLOV EXPONTNTIAL INTEGRATORS
   !
   !----------------------------------

   !run_test = .false.
   !if (run_test) then
   !   ntau = 20; allocate(tauv(1:ntau)); tauv = logspace(-5.0_wp,0.0_wp,ntau)
   !   torder = 1
   !   call run_kexpm_test(LTI%A, LTI%prop, U0(1), tauv, torder, 1000)
   !end if

   !----------------------------------
   !
   ! DLRA TEST FOR RICCATI EQUATION
   !
   !----------------------------------

   !run_test = .false.
   !if (run_test) then
   !   nrk  = 6; allocate(rkv(1:nrk));   rkv  = (/ 6, 10, 12, 14, 20, 40 /)
   !   ntau = 5; allocate(tauv(1:ntau)); tauv = (/ 1.0, 0.1, 0.01, 0.001, 0.0001 /)
   !   Tend = 1.0_wp
   !   nrep = 60
   !   ! run DLRA
   !   ifsave = .true. ! save X and Y matrices to disk (LightROM/local)
   !   ifverb = .true. ! verbosity
   !   iflogs = .true. ! write logs with convergence and signular value evolution
   !   call run_DLRA_riccati_test(LTI, U0, S0, Qc, Rinv, &
   !                              & rkv, tauv, Tend, nrep, & 
   !                              & ifsave, ifverb, iflogs)
   !   deallocate(rkv)
   !   deallocate(tauv)
   !end if 

   !----------------------------------
   !
   ! RANK-ADAPTIVE DLRA LYAPUNOV
   !
   !----------------------------------

   !run_test = .false.
   !if (run_test) then
   !   nrk  = 1; allocate(rkv(1:nrk));   rkv  = (/ 6 /)
   !   ntau = 1; allocate(tauv(1:ntau)); tauv = (/ 0.1 /)
   !   allocate(TOv(1)); TOv = (/ 1 /)
   !   Tend = 1.0_wp
   !   nrep = 60
   !   ! run DLRA
   !   ifsave = .false. ! save X and Y matrices to disk (LightROM/local)
   !   ifverb = .false. ! verbosity
   !   iflogs = .false. ! write logs with convergence and signular value evolution
   !   call run_DLRA_rank_adaptive_test(LTI, U0, S0, rkv, tauv, TOv, Tend, nrep, ifsave, ifverb, iflogs)
   !   deallocate(rkv)
   !   deallocate(tauv)
   !   deallocate(TOv)
   !end if 

   !----------------------------------
   !
   ! RANK-ADAPTIVE DLRA LYAPUNOV -- CONVERGENCE with increment norm
   !
   !----------------------------------

   !run_test = .true.
   !if (run_test) then
   !   nrk  = 1; allocate(rkv(1:nrk));   rkv  = (/ 6 /)
   !   ntau = 1; allocate(tauv(1:ntau)); tauv = (/ 0.1 /)
   !   allocate(TOv(1)); TOv = (/ 1 /)
   !   Tend = 150.0_wp
   !   ! run DLRA
   !   ifsave = .false. ! save X and Y matrices to disk (LightROM/local)
   !   ifverb = .false. ! verbosity
   !   iflogs = .false. ! write logs with convergence and signular value evolution
   !   call logger%configure(level=warning_level)
   !   call run_DLRA_rank_adaptive_test(LTI, U0, S0, rkv, tauv, TOv, Tend, 1, ifsave, ifverb, iflogs)
   !   deallocate(rkv)
   !   deallocate(tauv)
   !   deallocate(TOv)
   !end if 

   return
end program demo