program demo
   ! Standard Library
   use stdlib_optval, only : optval 
   use stdlib_linalg, only : eye, svdvals, eye
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

   character(len=128), parameter :: this_module = 'Laplacian2D_LTI_Lyapunov_Main'

   !----------------------------------------------------------
   !-----     LYAPUNOV EQUATION FOR LAPLACE OPERATOR     -----
   !----------------------------------------------------------

   ! DLRA
   integer, parameter :: rkmax = 16
   integer, parameter :: rk_X0 = 2
   logical, parameter :: verb  = .false.!.true.
   logical, parameter :: save  = .false.
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

   ! DLRA opts
   type(dlra_opts) :: opts



   ! CHK
   integer   :: is, ie, irow, istep, nstep, print_from, print_to
   logical   :: ifprint, ifprintsvd, iportho, ipdebug, ifpivot
   real(wp), allocatable  :: sv(:,:), R(:,:), wrk(:,:)
   integer, dimension(2) :: shape_out
   class(abstract_vector_rdp),  allocatable   :: exptAU
   class(abstract_vector_rdp),  allocatable   :: U1(:), Uwrk(:)
   class(abstract_vector_rdp),  allocatable   :: BBTU(:)
   character(len=128) :: pfx
   integer :: perm(rkmax)

   

   !call logger%configure(level=error_level, time_stamp=.false.); print *, 'Logging set to error_level.'
   call logger%configure(level=warning_level, time_stamp=.false.); print *, 'Logging set to warning_level.'

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
   Vdata = matmul(transpose(Q), matmul(-BBTdata, Q))
   Wdata = matmul(transpose(Z), matmul(   Vdata, Z))
   
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
   X0 = CALE(Xref, Adata, BBTdata)
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

  !   X = LR_state()
   Tend = 0.01_wp
   rk = 8
   dt = 0.002_wp
   nstep = nint(Tend/dt)
   print_from = 1
   print_to   = 2
   ifprint = .false.
   ifprintsvd = .true.
   iportho = .true.
   ipdebug = .true.
   ifpivot = .true.


   irow = 4
   shape_out = [ irow, int(rk/irow) ]
   
   call logger%configure(level=information_level, time_stamp=.false.); print *, 'Logging set to information_level.'
   call X%initialize_LR_state(U0(:rk_X0), S0(:rk_X0,:rk_X0), rk, rk, .false.)
   call logger%configure(level=error_level, time_stamp=.false.); print *, 'Logging set to error_level.'
   ! for safety
   call orthonormalize_basis(X%U)

   allocate(  U1(rk), source=X%U(1))
   allocate(BBTU(rk), source=X%U(1))

   call get_state(U_out(:,:rk), X%U(:rk), 'Get initial conditions')
   if (save) then
      call save_npy(trim(fldr)//'U0_full.npy', U_out(:,:rk))
      call save_npy(trim(fldr)//'S0_full.npy', X%S(:rk,:rk))
   end if
   print *, 'Read full IC from file:', 'U0_full.npy', 'S0_full.npy'
   call load_npy(trim(fldr)//'U0_full.npy', U_load)
   call set_state(X%U(:rk), U_load(:,:rk), 'Set initial conditions')
   call load_npy(trim(fldr)//'S0_full.npy', U_load)
   X%S(:rk,:rk) = U_load(:rk,:rk)

   ! print
   !call pvec(X%U, 'U0:')
   !call pmat(X%S, 'S0:')
   
   print *, ''
   print '(A)', 'Solver parameters'
   print '(3X,A,F6.4)', 'dt = ', dt
   print '(3X,A,I0)',   'TO = ', 1
   print '(3X,A,I0)',   'rk = ', rk
   print *, ''

   call pmat(reshape(svdvals(X%S(:rk,:rk)), shape_out), name='SVD(S0)', p=12, prefix='# INIT')
   allocate(exptAU, source=X%U(1)); call exptAU%zero()
   allocate(R(rk,rk));   R   = 0.0_wp
   allocate(wrk(rk,rk)); wrk = 0.0_wp

   do istep = 1, nstep
      if (istep .ge. print_from .and. istep .le. print_to) ifprint = .true.
      if (ifprintsvd .or. ipdebug .or. iportho) then
         print *, ''
         print *, ''
         print *, 'step ', istep
         print *, ''
         print *, ''
      !
      ! M step
      !
      !print *, ''
      !print *, '##'
      !print *, '##      M'
      !print *, '##'
      end if
      !print *, ''
      write(pfx,'(I4,3X,A)'), istep, ' M    '
      !
      !call M_forward_map_rdp(X, LTI%A, dt, info, exptA=exptA, iftrans=.false.)
      !
      if (ifprint .and. ipdebug) call pvec(X%U, 'U0:', prefix=pfx)

      do i = 1, rk
         call exptA(exptAU, LTI%A, X%U(i), dt, info, .false.)
         call X%U(i)%axpby(0.0_wp, exptAU, 1.0_wp) ! overwrite old solution
      end do

      if (ifprint .and. ipdebug) call pvec(X%U, 'U1:', prefix=pfx)

      ! Reorthonormalize in-place
      if (ifpivot) then
         do i = 1, rk
            print *, i, X%U(i)%dot(X%U(i))
         end do

         call qr(X%U(:rk), R, perm(:rk), info)
         call check_info(info, 'qr_pivot', module=this_module, procedure='M_forward_map_rdp')
         call apply_inverse_permutation_matrix(R, perm(:rk))
         print '(A,*(I0,1X))', 'permutation: ', perm(:rk)-1
      else
         call qr(X%U(:rk), R, info)
         call check_info(info, 'qr', module=this_module, procedure='M_forward_map_rdp')
      end if

      if (istep.eq.2) STOP 6

      if (ifprint .and. iportho) then
         call innerprod(wrk, X%U(:rk), X%U(:rk))
         print '(A,3X,A,E16.9)', trim(pfx), 'QR CHECK: Orthogonality        max(Q.T @ Q - I) =', maxval(wrk - eye(rk))
         call pvec(X%U(:rk), 'Q:', prefix=pfx)
         call pmat(R, 'R:', prefix=pfx)
         if (.not.ifpivot) then
            print *, 'Rii:'
            print '(A,3X,*(E16.9,1X))', trim(pfx), ( R(i,i), i=1,rk )
         end if
      end if

      ! Update coefficient matrix
      X%S(:rk,:rk) = matmul(R, matmul(X%S(:rk,:rk), transpose(R)))
      !
      if (ifprintsvd) then
         print *, ''
         call pmat(reshape(svdvals(X%S(:rk,:rk)), shape_out), 'SVD(X%S)', p=12, prefix=pfx)
      end if
      !
      ! G step
      !
      !print *, ''
      !print *, '##'
      !print *, '##      G'
      !print *, '##'
      !print *, ''
      !
      ! call G_forward_map_lyapunov(X, LTI%B, dt, info)
      !
      !print *, ''
      !print *, '       ##'
      !print *, '       ##      K'
      !print *, '       ##'
      !print *, ''
      write(pfx,'(I4,3X,A)'), istep, ' G    K    '
      call zero_basis(U1)
      call zero_basis(BBTU)
      call linear_combination(Uwrk, X%U(:rk), X%S(:rk,:rk))  ! K0
      call copy(U1, Uwrk)
      
      if (ifprint .and. ipdebug) call pvec(U1, 'K0:', prefix=pfx)

      call apply_outerprod(BBTU, LTI%B, X%U(:rk))                ! Kdot
      !call pvec(BBTU, 'Kdot:', prefix=pfx)
      ! Construct intermediate solution U1
      call axpby_basis(U1, 1.0_wp, BBTU, dt)           ! K0 + tau*Kdot
      !call pvec(U1, 'K1:', prefix=pfx)
      ! Orthonormalize in-place
      if (ifpivot) then
         call qr(U1, X%S(:rk,:rk), perm(:rk), info)
         call check_info(info, 'qr_pivot', module=this_module, procedure='M_forward_map_rdp')
         call apply_inverse_permutation_matrix(X%S(:rk,:rk), perm(:rk))
         print '(A,*(I0,1X))', 'permutation: ', perm(:rk)-1
      else
         call qr(U1, X%S(:rk,:rk), info, tol=atol_dp)
         call check_info(info, 'qr', module=this_module, procedure='K_step_Lyapunov_rdp')
      end if

      if (ifprint .and. iportho) then
         call innerprod(wrk, X%U(:rk), X%U(:rk))
         print '(A,3X,A,E16.9)', trim(pfx), 'QR CHECK: Orthogonality        max(Q.T @ Q - I) =', maxval(wrk - eye(rk))
         call pvec(U1, 'Q:', prefix=pfx)
         call pmat(X%S(:rk,:rk), 'Sh = R:', prefix=pfx)
         if (.not.ifpivot) then
            print *, 'Shii:'
            print '(A,3X,*(E16.9,1X))', trim(pfx), ( X%S(i,i), i=1,rk )
         end if
      end if

      if (ifprintsvd) then
         print *, ''
         call pmat(reshape(svdvals(X%S(:rk,:rk)), shape_out), 'SVD(X%S)', p=12, prefix=pfx)
      end if
      !print *, ''
      !print *, '       ##'
      !print *, '       ##      S'
      !print *, '       ##'
      !print *, ''
      write(pfx,'(I4,3X,A)'), istep, ' G    S    '
      call innerprod(wrk, U1, BBTU)          ! - Sdot
      ! Construct intermediate coefficient matrix
      X%S(:rk,:rk) = X%S(:rk,:rk) - dt*wrk

      if (ifprint .and. ipdebug) call pmat(X%S(:rk,:rk), 'S1:', prefix=pfx)

      if (ifprintsvd) then
         print *, ''
         call pmat(reshape(svdvals(X%S(:rk,:rk)), shape_out), 'SVD(X%S)', p=12, prefix=pfx)
      end if
      !print *, ''
      !print *, '       ##'
      !print *, '       ##      L'
      !print *, '       ##'
      !print *, ''
      write(pfx,'(I4,3X,A)'), istep, ' G    L    '
      
      if (ifprint .and. ipdebug) call pvec(X%U, 'UA:', prefix=pfx)

      !call pmat(transpose(X%S), 'St.T:', prefix=pfx)
      call linear_combination(Uwrk, X%U(:rk), transpose(X%S(:rk,:rk)))  ! L0.T
      !call pvec(Uwrk, 'L0.T:', prefix=pfx)
      ! Construct derivative
      call apply_outerprod(X%U(:rk), LTI%B, U1)       ! Ldot.T
      !call pvec(Uwrk, 'Ldot.T = B @ B.T @ U1:', prefix=pfx)
      ! Construct solution L1.T
      call axpby_basis(Uwrk, 1.0_wp, X%U(:rk), dt)
      !call pvec(Uwrk, 'L1.T = L0.T + dt*Ldot.T:', prefix=pfx)
      ! Update coefficient matrix
      call innerprod(X%S(:rk,:rk), Uwrk, U1)
      !call pmat(X%S(:rk,:rk), 'S:', prefix=pfx)

      call copy(X%U(:rk), U1)
      
      !print *, ''
      write(pfx,'(I4,3X,A)'), istep, ' G    '
      if (ifprintsvd) call pmat(reshape(svdvals(X%S(:rk,:rk)), shape_out), 'SVD(X%S)', p=12, prefix=pfx)
      ifprint = .false.
   end do
   print *, ''
   print *, ' ran integrator for ', nstep, 'steps.' 
   print *, ''
   call pmat(reshape(svdvals(X%S(:rk,:rk)), shape_out), name='SVD(S)', p=12, prefix='# EXIT')
   print *, ''
   print *, ''
   print *, ''
   print *, ''
   print *, ''
   print *, ''

contains
   subroutine pmat(A, name, row, col, p, prefix)
      real(wp), intent(in) :: A(:,:)
      character(len=*), optional, intent(in) :: name
      integer, optional, intent(in) :: row
      integer, optional, intent(in) :: col
      integer, optional, intent(in) :: p
      character(len=*), optional, intent(in) :: prefix
      ! internal
      integer :: i, row_, col_, p_, w
      character(len=128) :: name_, prefix_
      character(len=128) :: fmt

      row_ = optval(row, size(A,1))
      col_ = optval(col, size(A,2))
      p_   = optval(p, 5)
      prefix_ = optval(prefix, '')
      w = p_ + 3
      write(fmt,'("(A,3X,*(F",I0,".",I0,",1X))")') w, p_
      if (row_ .gt. size(A,1)) call logger%log_error('row > row(A)')
      if (col_ .gt. size(A,2)) call logger%log_error('col > col(A)')
      name_ = optval(name, 'mat')

      print *, name_
      do i = 1,row_
         print fmt, trim(prefix_), A(i,:col_)
      end do
      return
   end subroutine pmat

   subroutine pvec(U, name, col, p, prefix)
      class(abstract_vector_rdp), intent(in) :: U(:)
      character(len=*), optional, intent(in) :: name
      integer, optional, intent(in) :: col
      integer, optional, intent(in) :: p
      character(len=*), optional, intent(in) :: prefix
      ! internal
      real(wp), allocatable :: Umat(:,:)
      integer :: i, col_, kdim
      character(len=128) :: name_

      select type (U)
      type is (state_vector)
         kdim = size(U)
         col_ = optval(col, kdim)
         if (col_ .gt. kdim) call logger%log_error('col > col(U)')
         name_ = optval(name, 'mat')
         allocate(Umat(N,col_))
         call get_state(Umat, U(:col_), name)
         call pmat(Umat, name=name_, row=N, col=col_, p=p, prefix=prefix)
      end select
      return     
   end subroutine pvec
end program demo