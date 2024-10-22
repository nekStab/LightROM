program demo
   ! Standard Library
   use stdlib_optval, only : optval 
   use stdlib_linalg, only : eye, svdvals
   use stdlib_math, only : all_close, logspace
   use stdlib_io_npy, only : save_npy, load_npy
   use stdlib_logger, only : information_level, warning_level, debug_level, error_level, none_level
   ! LightKrylov for linear algebra
   use LightKrylov
   use LightKrylov, only : wp => dp
   use LightKrylov_Logger
   use LightKrylov_ExpmLib
   use LightKrylov_Utils
   ! LightROM
   use LightROM_AbstractLTIsystems
   use LightROM_Utils
   use LightROM_LyapunovSolvers
   use LightROM_LyapunovUtils
   ! Laplacian
   use Laplacian2D_LTI_Lyapunov_Base
   use Laplacian2D_LTI_Lyapunov_Operators
   use Laplacian2D_LTI_Lyapunov_RKlib
   use Laplacian2D_LTI_Lyapunov_Utils
   implicit none

   character*128, parameter :: this_module = 'Laplacian2D_LTI_Lyapunov_Main'

   !----------------------------------------------------------
   !-----     LYAPUNOV EQUATION FOR LAPLACE OPERATOR     -----
   !----------------------------------------------------------

   ! DLRA
   integer, parameter :: rkmax = 16
   integer, parameter :: rk_X0 = 2
   logical, parameter :: verb  = .false.!.true.
   logical, parameter :: save  = .false.!.true.
   logical, parameter :: read  = .false.
   character*128      :: oname, fldr
   ! rk_B is set in laplacian2D.f90

   integer  :: nrk, ndt, rk,  torder
   real(wp) :: dt, Tend
   ! vector of dt values
   real(wp), allocatable :: dtv(:)
   ! vector of rank values
   integer,  allocatable :: rkv(:)

   ! Exponential propagator (RKlib).
   type(rklib_lyapunov_mat), allocatable :: RK_propagator

   ! LTI system
   type(lti_system)                :: LTI
   real(wp), allocatable           :: D(:,:)
   integer                         :: p

   ! Laplacian
   type(laplace_operator), allocatable :: A

   ! LR representation
   type(LR_state)                  :: X
   type(state_vector), allocatable :: U(:)
   real(wp) , allocatable          :: S(:,:)
   
   !> STATE MATRIX (RKlib)
   type(state_matrix)              :: X_mat_RKlib(2)
   real(wp), allocatable           :: X_RKlib(:,:,:)
   real(wp), allocatable           :: X_DLRA(:,:,:)
   real(wp)                        :: X_RKlib_ref(N,N)

   ! Initial condition
   type(state_vector)              :: U0(rkmax)
   real(wp)                        :: S0(rkmax,rkmax)
   ! Matrix
   real(wp)                        :: U0_in(N,rkmax)
   real(wp)                        :: X0(N,N)

   ! OUTPUT
   real(wp)                        :: U_out(N,rkmax)
   real(wp)                        :: X_out(N,N)
   real(wp), allocatable           :: U_load(:,:)

   !> Information flag.
   integer                         :: info
   integer                         :: i, j, k, irep, nrep

   ! PROBLEM DEFINITION
   real(wp)  :: Adata(N,N)
   real(wp)  :: Bdata(N,N)
   real(wp)  :: BBTdata(N,N)
   
   ! LAPACK SOLUTION
   real(wp)  :: Xref(N,N)
   ! DSYTD2
   real(wp)  :: Dm(N), work(N), wr(N), wi(N)
   real(wp)  :: E(N-1), tw(N-1)
   real(wp)  :: T(N,N), Q(N,N), Z(N,N), Vdata(N,N), Wdata(N,N), Ydata(N,N)
   real(wp)  :: scale
   real(wp)  :: svals(N)
   integer   :: isgn

   ! timer
   integer   :: clock_rate, clock_start, clock_stop
   integer   :: is, ie, irow

   ! DLRA opts
   type(dlra_opts) :: opts

   !call logger%configure(level=error_level, time_stamp=.false.); print *, 'Logging set to error_level.'
   call logger%configure(level=warning_level, time_stamp=.false.); print *, 'Logging set to error_level.'

   call system_clock(count_rate=clock_rate)

   irow = 8

   print *, '---------------------------------------------'
   print *, '   DYNAMIC LOW-RANK APPROXIMATION  -  DLRA'
   print *, '---------------------------------------------'
   print *, ''
   print *, 'LYAPUNOV EQUATION FOR THE 2D LAPLACE OPERATOR:'
   print *, ''
   print *, '    Algebraic Lyapunov equation:'
   print *, '                0 = A @ X + X @ A.T + Q'
   print *, ''
   print *, '    Differential Lyapunov equation:'
   print *, '          \dot{X} = A @ X + X @ A.T + Q'
   print *, ''
   write(*,'(A16,I4,"x",I4)') '  Problem size: ', N, N
   print *, ''
   print *, '  Initial condition: rank(X0) =', rk_X0
   print *, '  Inhomogeneity:     rank(Q)  =', rk_B
   print *, ''
   print *, '---------------------------------------------'
   print *, ''

   fldr = 'example/DLRA_laplacian2D_lti_lyapunov/'

   ! Define RHS B
   do i = 1, rk_b
      call B(i)%rand(ifnorm = .false.)   
   end do
   call get_state(Bdata(:,:rk_b), B, 'Get B')
   if (save) call save_npy(trim(fldr)//'B.npy', Bdata)
   if (read) then
      call load_npy(trim(fldr)//'B.npy', U_load)
      call set_state(B, U_load(:,:rk_b), 'Load B')
      call get_state(Bdata(:,:rk_b), B, 'Get B')
   end if
   BBTdata = matmul(Bdata(:,:rk_b), transpose(Bdata(:,:rk_b)))
   BBT(:N**2) = reshape(BBTdata, shape(BBT))

   ! Define LTI system
   LTI = lti_system()
   call LTI%initialize_lti_system(A, B, B)
   call zero_basis(LTI%CT)

   ! Define initial condition of the form X0 + U0 @ S0 @ U0.T SPD
   if (verb) print *, '    Define initial condition'
   call generate_random_initial_condition(U0(:rk_X0), S0(:rk_X0,:rk_X0), rk_X0)
   call get_state(U_out(:,:rk_X0), U0(:rk_X0), 'Get initial conditions')
   if (save) then
      call save_npy(trim(fldr)//'U0.npy', U_out(:,:rk_X0))
      call save_npy(trim(fldr)//'S0.npy', S0)
   end if
   if (read) then
      call load_npy(trim(fldr)//'U0.npy', U_load)
      call zero_basis(U0); S0 = 0.0_wp
      call set_state(U0(:rk_X0), U_load(:,:rk_X0), 'Load U0')
      call load_npy(trim(fldr)//'S0.npy', U_load)
      S0(:rk_X0,:rk_X0) = U_load(:rk_X0,:rk_X0)
   end if
   
   ! Compute the full initial condition X0 = U_in @ S0 @ U_in.T
   X0 = matmul( U_out(:,:rk_X0), matmul(S0(:rk_X0,:rk_X0), transpose(U_out(:,:rk_X0))))
   if (save) call save_npy(trim(fldr)//'X0.npy', X0)
   if (read) then
      call load_npy(trim(fldr)//'X0.npy', U_load)
      X0 = U_load
   end if
   
   print *, 'SVD X0'
   svals = svdvals(X0)
   do i = 1, ceiling(rkmax*1.0/irow)
      is = (i-1)*irow+1; ie = min(i*irow, rkmax)
      print '(2X,I2,A,I2,*(F16.12,1X))', is, '-', ie, svals(is:ie)
   end do
   print *, ''
   
   !------------------
   ! COMPUTE EXACT SOLUTION OF THE LYAPUNOV EQUATION WITH LAPACK
   !------------------

   print *, 'I.   Exact solution of the algebraic Lyapunov equation (LAPACK):'
   call system_clock(count=clock_start)     ! Start Timer

   ! Explicit 2D laplacian
   call build_operator(Adata)
   if (save) call save_npy(trim(fldr)//'A.npy', Adata)
   if (read) then
      call load_npy(trim(fldr)//'A.npy', U_load)
      Adata = U_load
   end if

   ! Transform operator to tridiagonal form
   call dsytd2('L', N, Adata, N, Dm, E, tw, info)
   
   ! Reconstruct T and Q
   call reconstruct_TQ(T, Q, Adata(:N,:N), Dm, E, tw)
   
   ! compute real Schur form A = Z @ T @ Z.T
   call dhseqr('S', 'I', N, 1, N, T, N, wr, wi, Z, N, work, N, info )

   ! Change RHS Basis: base --> Q --> Z
   Vdata = matmul(transpose(Q), matmul(BBTdata, Q))
   Wdata = matmul(transpose(Z), matmul(  Vdata, Z))
   
   ! Compute solution of Lyapunov equation for Schur decomposition
   isgn = 1; scale = 0.1_wp
   call dtrsyl('N', 'T', isgn, N, N, T, N, T, N, Wdata, N, scale, info)

   ! Return to original basis to obtain X_ref: Z --> Q --> base
   Ydata = matmul(Z, matmul(Wdata, transpose(Z)))
   Xref  = matmul(Q, matmul(Ydata, transpose(Q)))

   call system_clock(count=clock_stop)      ! Stop Timer
   write(*,'(A40,F10.4," s")') '--> X_ref.    Elapsed time:', real(clock_stop-clock_start)/real(clock_rate)
   print *, ''
   if (save) call save_npy(trim(fldr)//'Xref.npy', Xref)
   if (read) then
      call load_npy(trim(fldr)//'Xref.npy', U_load)
      Xref = U_load
   end if

   call build_operator(Adata)
   ! sanity check
   X0 = CALE(Xref, Adata, -BBTdata)
   print *, '    Direct problem: || res(X_ref) ||_2/N', norm2(X0)/N
   print *, ''
   ! compute svd
   print *, 'SVD Xref'
   svals = svdvals(Xref)
   do i = 1, ceiling(N*1.0_wp/irow)
      is = (i-1)*irow+1; ie = min(i*irow, N)
      print '(2X,I2,A,I2,*(F16.12,1X))', is, '-', ie, svals(is:ie)
   end do
   print *, ''

   !------------------
   ! COMPUTE SOLUTION WITH RK FOR DIFFERENT INTEGRATION TIMES AND COMPARE TO STUART-BARTELS
   !------------------

   print *, 'II.  Compute solution of the differential Lyapunov equation (Runge-Kutta):'
   print *, ''
   ! initialize exponential propagator
   nrep = 10
   Tend = 0.001_wp
   RK_propagator = rklib_lyapunov_mat(Tend)

   allocate(X_RKlib(N, N, nrep))
   call get_state(U_out(:,:rk_X0), U0(:rk_X0), 'Get initial conditions')
   X0 = matmul( U_out(:,:rk_X0), matmul(S0(:rk_X0,:rk_X0), transpose(U_out(:,:rk_X0))))
   call set_state(X_mat_RKlib(1:1), X0, 'Set RK X0')
   write(*,'(A10,A26,A26,A20)') 'RKlib:','Tend','|| X_RK - X_ref ||_2/N', 'Elapsed time'
   print *, '         ------------------------------------------------------------------------'
   do irep = 1, nrep
      call system_clock(count=clock_start)     ! Start Timer
      ! integrate
      call RK_propagator%matvec(X_mat_RKlib(1), X_mat_RKlib(2))
      ! recover output
      call get_state(X_RKlib(:,:,irep), X_mat_RKlib(2:2), 'Get RK solution')
      ! replace input
      call set_state(X_mat_RKlib(1:1), X_RKlib(:,:,irep), 'Reset RK X0')
      call system_clock(count=clock_stop)      ! Stop Timer
      write(*,'(I10,F26.4,E26.8,F18.4," s")') irep, irep*Tend, &
                     & norm2(X_RKlib(:,:,irep) - Xref)/N, &
                     & real(clock_stop-clock_start)/real(clock_rate)
   end do

   ! Choose relevant reference case from RKlib
   X_RKlib_ref = X_RKlib(:,:,1)
   if (save) call save_npy(trim(fldr)//'Xref_RK.npy', X_RKlib_ref)
   if (read) then
      call load_npy(trim(fldr)//'Xref_RK.npy', U_load)
      X_RKlib_ref = U_load
   end if

   print *, ''
   print *, 'SVD Xrk'
   svals = svdvals(X_RKlib_ref)
   do i = 1, ceiling(N*1.0_wp/irow)
      is = (i-1)*irow+1; ie = min(i*irow, N)
      print '(2X,I2,A,I2,*(F16.12,1X))', is, '-', ie, svals(is:ie)
   end do
   print *, ''
   
   !------------------
   ! COMPUTE DLRA FOR SHORTEST INTEGRATION TIMES FOR DIFFERENT DT AND COMPARE WITH RK SOLUTION
   !------------------

   print *, ''
   print *, 'III. Compute approximate solution of the differential Lyapunov equation using DLRA:'
   print *, ''
   write(*,'(A10,A4,A4,A10,A8,A26,A26,A20)') 'DLRA:','  rk',' TO','dt','Tend','|| X_DLRA - X_RK ||_2/N','|| X_DLRA - Xref ||_2/N', 'Elapsed time'
   print *, '         ------------------------------------------------------------------------'









































   
   ! Choose input ranks and integration steps
   rkv = [ 2, 6, 10, 14, 16 ]
   dtv = logspace(-6.0_wp, -3.0_wp, 4, 10)
   dtv = dtv(size(dtv):1:-1)

   allocate(X_DLRA(N, N, size(dtv)*size(rkv)*2))

   irep = 0
   X = LR_state()
   do i = 1, size(rkv)
      rk = rkv(i)
      do torder = 1, 2
         do j = 1, size(dtv)
            irep = irep + 1
            dt = dtv(j)

            ! Reset input
            call X%initialize_LR_state(U0, S0, rk)

            ! run step
            opts = dlra_opts(mode=torder, if_rank_adaptive=.false.)
            call system_clock(count=clock_start)     ! Start Timer
            call projector_splitting_DLRA_lyapunov_integrator(X, LTI%A, LTI%B, Tend, dt, info, exptA=exptA, options=opts)
            call system_clock(count=clock_stop)      ! Stop Timer

            ! Reconstruct solution
            call get_state(U_out(:,:rk), X%U, 'Reconstruct solution')
            X_out = matmul(U_out(:,:rk), matmul(X%S, transpose(U_out(:,:rk))))
            write(oname,'("data_X_DRLA_TO",I1,"_rk",I2.2,"_",I3.3,".npy")') torder, rk, irep
            if (save) call save_npy(trim(fldr)//oname, X_out)
            if (read) then
               call load_npy(trim(fldr)//oname, U_load)
               X_out = U_load
            end if
      
            write(*,'(A10,I4," TO",I1,F10.6,F8.4,E26.8,E26.8,F18.4," s")') 'OUTPUT', &
                              & rk, torder, dt, Tend, &
                              & norm2(X_RKlib_ref - X_out)/N, norm2(X_out - Xref)/N, &
                              & real(clock_stop-clock_start)/real(clock_rate)

            deallocate(X%U)
            deallocate(X%S)
            X_DLRA(:,:,irep) = X_out
         end do
         print *, ''
      end do
      print *, ''
   end do
   nrep = irep

   print *, 'SVD X_DLRA final'
   svals = svdvals(X_out)
   do i = 1, ceiling(N*1.0_wp/irow)
      is = (i-1)*irow+1; ie = min(i*irow, N)
      print '(2X,I2,A,I2,*(F16.12,1X))', is, '-', ie, svals(is:ie)
   end do
   print *, ''

   print *, 'nrep', nrep, 'X_DLRA', size(X_DLRA)

   if (save) then
      call save_npy("example/DLRA_laplacian2D/data_A.npy", Adata)
      call save_npy("example/DLRA_laplacian2D/data_X0.npy", X0)
      call save_npy("example/DLRA_laplacian2D/data_Q.npy", BBTdata)
      call save_npy("example/DLRA_laplacian2D/data_X_ref.npy", Xref)
      call save_npy("example/DLRA_laplacian2D/data_X_RK.npy", X_RKlib_ref)
   end if

end program demo