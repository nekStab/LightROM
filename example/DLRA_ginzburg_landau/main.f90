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

   character*128      :: onameU, onameS, oname
   ! rk_B & rk_C are set in ginzburg_landau_base.f90

   integer, parameter :: irow = 8
   integer  :: nrk, ntau, rk,  torder, nsteps, iref
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
   type(dlra_opts)                           :: opts

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
   real(wp)                                  :: res(N,N)

   ! Reference solutions (BS & RK)
   real(wp)                                  :: Xref_BS(N,N)
   real(wp)                                  :: Xref_RK(N,N)

   ! IO
   real(wp),           allocatable           :: U_load(:,:)
   real(wp),           allocatable           :: BBTW_load(:,:)

   ! Information flag.
   integer                                   :: info

   ! Counters
   integer                                   :: i, j, k, irep, nrep, istep, is, ie
   integer,                allocatable       :: perm(:)
   real(wp)                                  :: etime

   real(wp)                           :: U_svd(2*nx,2*nx)
   real(wp)                           :: S_svd(2*nx)
   real(wp)                           :: V_svd(2*nx,2*nx)
   real(wp), dimension(:,:),                   allocatable :: svals(:)

   real(wp) :: wrk(2,1)

   logical        :: ifsave, ifload, iflogs, ifW, ifchk
   
   ! Exact comparison
   type(LR_state),                allocatable   :: X
   class(abstract_vector_rdp),             allocatable   :: exptAU    ! scratch basis
   real(wp),                               allocatable   :: R(:,:), S0ref(:,:)    ! QR coefficient matrix
   class(abstract_vector_rdp),  allocatable   :: U0ref(:)
   class(abstract_vector_rdp),  allocatable   :: BBTU(:)
   real(wp), allocatable :: ssvd_v(:)
   character(len=128) :: fmt
   real(wp), allocatable :: xv(:)
   real(wp), dimension(:,:),                   allocatable :: ssvd_r
   integer, dimension(2) :: shape_out
   integer :: istart, iend, ndt, iostat
   character(len=128), allocatable :: iomsg

   !call logger%configure(level=information_level, time_stamp=.false.)
   call logger%configure(level=error_level, time_stamp=.false.)

   !----------------------------------
   !-----     INITIALIZATION     -----
   !----------------------------------

   ! Initialize mesh and system parameters A, B, CT
   print *, 'Initialize parameters'
   call initialize_parameters()

   ! Initialize propagator
   print *, 'Initialize exponential propagator'
   prop = exponential_prop(1.0_wp)

   ! Initialize LTI system
   A = GL_operator()
   print *, 'Initialize LTI system (A, prop, B, CT, _)'
   LTI = lti_system()
   call LTI%initialize_lti_system(A, prop, B, CT)

   svals = svdvals(BBTW)
   print *, ''
   print '(1X,A,*(F16.12,X))', 'SVD(1:3) BBTW:   ', svals(1:3)

   print *, ''
   print *, 'Check residual computation with Bartels-Stuart solution:'
   oname = trim(basepath)//"BS/CGL_Lyapunov_Controllability_Xref_BS_W.npy"
   call load_data(oname, U_load)
   Xref_BS = U_load
   
   print *, ''
   print '(A,F16.12)', '  |  X_BS  |/N = ', norm2(Xref_BS)/N
   print '(A,F16.12)', '  | res_BS |/N = ', norm2(CALE(Xref_BS, BBTW, .false.))/N
   print *, ''
   ! compute svd
   !call svd(Xref_BS, S_svd, U_svd, V_svd)
   svals = svdvals(Xref_BS)
   print *, 'SVD Xref:'
   do i = 1, ceiling(60.0/irow)
      is = (i-1)*irow+1; ie = i*irow
      print '(2X,I2,A,I2,*(1X,F16.12))', is, '-', ie, ( svals(j), j = is, ie )
   end do
   print *, ''
   
   ! Define initial condition of the form X0 = U0 @ S0 @ U0.T (SPD) or load from file
   allocate(U0(rk_X0), source=B(1)); call zero_basis(U0)
   allocate(S0(rk_X0,rk_X0)); S0 = 0.0_wp
   print *, 'Define initial condition'
   call generate_random_initial_condition(U0, S0, rk_X0)
   call get_state(U_out(:,:rk_X0), U0)
   X_out = matmul(matmul( U_out(:,:rk_X0), matmul( S0(:rk_X0,:rk_X0), transpose(U_out(:,:rk_X0)) ) ), weight_mat)
   print *, ''
   print '(A,F16.12)', '  |  X_0  |/N = ', norm2(X_out)/N
   print '(A,F16.12)', '  | res_0 |/N = ', norm2(CALE(X_out, BBTW, .false.))/N
   print *, ''
   ! compute svd
   svals = svdvals(X_out)
   do i = 1, ceiling(rk_X0*1.0_wp/irow)
      is = (i-1)*irow+1; ie = i*irow
      print '(2X,I2,A,I2,*(F16.12,1X))', is, '-', ie, ( svals(j), j = is, ie )
   end do
   print *, ''
   
   ! Run RK integrator for the Lyapunov equation
   T_RK   = 0.01_wp
   nsteps = 10
   iref   = 1
   call run_lyap_reference_RK(LTI, Xref_BS, Xref_RK, U0, S0, T_RK, nsteps, iref)

   print *, ''
   print '(A,F16.12)', '  |  X_RK  |/N = ', norm2(Xref_BS)/N
   print '(A,F16.12)', '  | res_RK |/N = ', norm2(CALE(Xref_BS, BBTW, .false.))/N
   print *, ''
   
   ! basic settings
   Tend = T_RK/nsteps*iref
   ifsave = .false. ! save data to disk (LightROM/local)

   ! DLRA with fixed rank
   rkv = [ 4, 10, 20 ]
   tauv = logspace(-5.0_wp, -3.0_wp, 3, 10)
   tauv = tauv(size(tauv):1:-1) ! reverse vector
   TOv  = [ 1, 2 ] 
   call run_lyap_DLRA_test(LTI, Xref_BS, Xref_RK, U0, S0, Tend, tauv, rkv, TOv, ifsave)
   
   ! DLRA with adaptive rank
   !tauv = logspace(-5.0_wp, -2.0_wp, 4, 10)
   !tauv = tauv(size(tauv):1:-1) ! reverse vector
   !TOv  = [ 2 ] ![ 1, 2 ]
   !tolv = [ 1e-10 ] ! [ 1e-2_wp, 1e-6_wp, 1e-10_wp ]
   !call run_lyap_DLRArk_test(LTI, Xref_BS, Xref_RK, U0, S0, Tend, tauv, TOv, tolv, ifsave)

   return
end program demo