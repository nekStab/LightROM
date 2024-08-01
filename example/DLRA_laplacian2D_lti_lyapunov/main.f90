program demo
   ! Standard Library
   use stdlib_optval, only : optval 
   use stdlib_linalg, only : eye, svdvals
   use stdlib_math, only : all_close, logspace
   use stdlib_io_npy, only : save_npy
   use stdlib_logger, only : error_level, none_level
   ! LightKrylov for linear algebra
   use LightKrylov
   use LightKrylov, only : wp => dp
   use LightKrylov_Constants
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
   integer, parameter :: rkmax = 20
   integer, parameter :: rk_X0 = 10
   logical            :: verb      = .false.
   logical            :: if_save   = .false. ! save outputs
   logical            :: if_convRK = .false. ! run RK to convergence or just enough for the test
   logical            :: if_rkad   = .false. 
   logical            :: if_kexpm  = .false.
   logical            :: if_t_all  = .true.
   character*128      :: oname
   ! rk_B is set in laplacian2D.f90

   integer  :: ndt, rk,  torder
   real(wp) :: dt, Tend, tchk
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
   type(rklib_exptA_laplacian), allocatable :: prop

   ! LR representation
   type(LR_state)                  :: X
   type(state_vector), allocatable :: U(:)
   real(wp) , allocatable          :: S(:,:)
   
   !> STATE MATRIX (RKlib)
   type(state_matrix)              :: X_mat_RKlib(2)
   real(wp), allocatable           :: X_RKlib(:,:,:)
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

   !> Information flag.
   integer                         :: info
   integer                         :: i, j, k, irep, nrep, rki

   ! PROBLEM DEFINITION
   real(wp)  :: Adata(N,N)
   real(wp)  :: Bdata(N,rkmax)
   real(wp)  :: BBTdata(N,N)
   
   ! LAPACK SOLUTION
   ! DSYTD2
   real(wp), dimension(N)   :: Dm, work, wr, wi
   real(wp), dimension(N-1) :: E, tw
   real(wp), dimension(N,N) :: Xref, T, Q, Z, Vdata, Wdata, Ydata
   real(wp)  :: scale
   integer   :: isgn
   ! SVD
   real(wp), allocatable :: svals(:)

   ! timer
   integer   :: clock_rate, clock_start, clock_stop
   real(wp)  :: etime

   ! DLRA opts
   type(dlra_opts) :: opts

   !----------------------------------------------------
   ! TEST PARAMETERS
   !----------------------------------------------------
   !
   ! Rank of initial condition
   rkv = (/ 14 /)
   ! Time Step
   ndt = 5
   dtv = logspace(-5.0_wp, -1.0_wp, ndt)
   ! Flags
   verb      = .false. ! verbosity
   if_save   = .false. ! save outputs
   if_convRK = .false. ! run RK to convergence or just enough for the test
   if_rkad   = .false. ! run rank adaptive
   if_kexpm  = .true.  ! use kexpm for the action of the matrix exponential (or RKlib)
   ! override
   if_t_all  = .true.  ! test all 4 options (+- rank adaptive, kexpm + RKlib)
   !----------------------------------------------------

   call logger_setup()
   call logger%configure(level=error_level, time_stamp=.false.); write(*,*) 'Logging set to error_level.'
   call system_clock(count_rate=clock_rate)

   write(*,*) '---------------------------------------------'
   write(*,*) '   DYNAMIC LOW-RANK APPROXIMATION  -  DLRA'
   write(*,*) '---------------------------------------------'
   write(*,*)
   write(*,*) 'LYAPUNOV EQUATION FOR THE 2D LAPLACE OPERATOR:'
   write(*,*)
   write(*,*) '    Algebraic Lyapunov equation:'
   write(*,*) '                0 = A @ X + X @ A.T + Q'
   write(*,*)
   write(*,*) '    Differential Lyapunov equation:'
   write(*,*) '          \dot{X} = A @ X + X @ A.T + Q'
   write(*,*)
   write(*,'(A16,I4,"x",I4)') '  Problem size: ', N, N
   write(*,*)
   write(*,*) '  Initial condition: rank(X0) =', rk_X0
   write(*,*) '  Inhomogeneity:     rank(Q)  =', rk_B
   write(*,*)
   write(*,*) '---------------------------------------------'
   write(*,*)

   ! Define RHS B
   call init_rand(B, ifnorm = .false.)
   call get_state(Bdata(:,:rk_b), B)
   BBTdata = -matmul(Bdata(:,:rk_b), transpose(Bdata(:,:rk_b)))
   BBT(:N**2) = -reshape(BBTdata, shape(BBT))

   ! Define LTI system
   prop = rklib_exptA_laplacian(1.0_wp)
   call LTI%initialize(A, prop, B, B)
   call zero_basis(LTI%CT)

   ! Define initial condition of the form X0 + U0 @ S0 @ U0.T SPD 
   if (verb .and. io_rank()) write(*,*) '    Define initial condition'
   call generate_random_initial_condition(U0(:rk_X0), S0(:rk_X0,:rk_X0), rk_X0)
   call get_state(U_out(:,:rk_X0), U0(:rk_X0))
   
   ! Compute the full initial condition X0 = U_in @ S0 @ U_in.T
   X0 = matmul( U_out(:,:rk_X0), matmul(S0(:rk_X0,:rk_X0), transpose(U_out(:,:rk_X0))))
   
   !------------------
   ! COMPUTE EXACT SOLUTION OF THE LYAPUNOV EQUATION WITH LAPACK
   !------------------

   if (io_rank()) write(*,*) 'I.   Exact solution of the algebraic Lyapunov equation (LAPACK):'
   call system_clock(count=clock_start)     ! Start Timer

   ! Explicit 2D laplacian
   call build_operator(Adata)

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
   etime = real(clock_stop-clock_start)/real(clock_rate)
   if (io_rank()) write(*,'(A40,F10.4," s")') '--> X_ref.    Elapsed time:', etime
   if (io_rank()) write(*,*)

   ! sanity check
   call build_operator(Adata) ! rebuild overwritten operator
   if (io_rank()) then
       write(*,*) '    Direct problem: || res ||_2/N =', norm2(CALE(Xref, Adata, BBT, .false.))/N
       write(*,*)
       call print_svdvals(Xref, 'Xref')
       write(*,*)
   end if
   
   !------------------
   ! COMPUTE SOLUTION WITH RK FOR DIFFERENT INTEGRATION TIMES AND COMPARE TO STUART-BARTELS
   !------------------

   if (io_rank()) write(*,*) 'II.  Compute solution of the differential Lyapunov equation (Runge-Kutta):'
   if (io_rank()) write(*,*)
   ! initialize exponential propagator
   Tend = 0.1_wp
   RK_propagator = rklib_lyapunov_mat(Tend)
   nrep = 1
   if (if_convRK) nrep = 10

   allocate(X_RKlib(N, N, nrep))
   ! set initial condition
   call set_state(X_mat_RKlib(1:1), X0)
   if (io_rank()) write(*,'(A17,A26,A26,A20)') 'RKlib:','Tend','|| X_RK - X_ref ||_2/N', 'Elapsed time'
   if (io_rank()) write(*,*) '         ------------------------------------------------------------------------'
   do irep = 1, nrep
      call system_clock(count=clock_start)     ! Start Timer
      ! integrate
      call RK_propagator%matvec(X_mat_RKlib(1), X_mat_RKlib(2))
      ! recover output
      call get_state(X_RKlib(:,:,irep), X_mat_RKlib(2:2))
      ! replace input
      call set_state(X_mat_RKlib(1:1), X_RKlib(:,:,irep))
      call system_clock(count=clock_stop)      ! Stop Timer
      etime = real(clock_stop-clock_start)/real(clock_rate)
      if (io_rank()) write(*,'("  OUTPUT_X   RK ",I2,F26.4,E26.8,F18.4," s")') irep, irep*Tend, &
                     & norm2(X_RKlib(:,:,irep) - Xref)/N, etime
      call print_svdvals(X_RKlib(:,:,irep), 'X_RK ')
      if (io_rank()) write(*,*)
   end do
 
   !------------------
   ! COMPUTE DLRA FOR SHORTEST INTEGRATION TIMES FOR DIFFERENT DT AND COMPARE WITH RK SOLUTION
   !------------------

   if (io_rank()) then
       write(*,*)
       write(*,*) 'III. Compute approximate solution of the differential Lyapunov equation using DLRA:'
       write(*,*)
       write(*,'(A10,A4,A4,A10,A8,A26,A20)') 'DLRA:','  rk',' TO','dt','Tend','|| X_DLRA - X_RK ||_2/N', 'Elapsed time'
       write(*,*) '       ------------------------------------------------------------------------'
   end if

   ! Choose relevant reference case from RKlib
   X_RKlib_ref = X_RKlib(:,:,1)

   do torder = 1, 2
      do i = 1, size(rkv)
         rk = rkv(i)
         if (verb .and. io_rank()) write(*,'(A10,I1)') ' torder = ', torder
         do j = ndt, 1, -1
            dt = dtv(j)
            if (verb .and. io_rank()) write(*,*) '    dt = ', dt, 'Tend = ', Tend
            tchk = 1.0_wp
            if (verb) tchk = max(0.01_wp, dt)
            ! run each option or only one of them
            if ((.not. if_rkad .and. .not. if_kexpm) .or. if_t_all) then
                if (io_rank()) call print_info(if_rkad, if_kexpm, Tend, dt, rk, torder)
                ! set DLRA opts
                opts = dlra_opts(mode=torder, verbose=verb, if_rank_adaptive=if_rkad, chktime=tchk)
                ! Initialize solution
                call X%initialize_LR_state(U0, S0, rk, rkmax)
                ! Run integrator
                call system_clock(count=clock_start)     ! Start Timer
                call projector_splitting_DLRA_lyapunov_integrator(X, LTI%prop, LTI%B, Tend, dt, info, & 
                       & exptA=exptA_rklib, options=opts)
                call system_clock(count=clock_stop)      ! Stop Timer
                etime = real(clock_stop-clock_start)/real(clock_rate)
                ! Reconstruct solution
                call get_state(U_out(:,:rk), X%U(:rk))
                X_out = matmul(U_out(:,:rk), matmul(X%S(:rk,:rk), transpose(U_out(:,:rk))))
                ! Print info
                if (io_rank()) then
                    write(*,'(A17,I4," TO",I1,F10.6,F8.4,E20.8,E20.8,F18.4," s")') 'OUTPUT_X   DLRA', &
                        & rk, torder, dt, Tend, norm2(X_RKlib_ref - X_out)/N, norm2(Xref - X_out)/N, etime
                    call print_svdvals(X_out, 'DLRA')
                end if
            else if ((.not. if_rkad .and. if_kexpm) .or. if_t_all) then
                if (io_rank()) call print_info(if_rkad, if_kexpm, Tend, dt, rk, torder)
                ! set DLRA opts
                opts = dlra_opts(mode=torder, verbose=verb, if_rank_adaptive=if_rkad, chktime=tchk)
                ! Initialize solution
                call X%initialize_LR_state(U0, S0, rk, rkmax)
                ! Run integrator
                call system_clock(count=clock_start)     ! Start Timer
                call projector_splitting_DLRA_lyapunov_integrator(X, LTI%A,    LTI%B, Tend, dt, info, &
                       & exptA=exptA,       options=opts)
                call system_clock(count=clock_stop)      ! Stop Timer
                etime = real(clock_stop-clock_start)/real(clock_rate)
                ! Reconstruct solution
                call get_state(U_out(:,:rk), X%U(:rk))
                X_out = matmul(U_out(:,:rk), matmul(X%S(:rk,:rk), transpose(U_out(:,:rk))))
                ! Print info
                if (io_rank()) then
                    write(*,'(A17,I4," TO",I1,F10.6,F8.4,E20.8,E20.8,F18.4," s")') 'OUTPUT_X   DLRA', &
                        & rk, torder, dt, Tend, norm2(X_RKlib_ref - X_out)/N, norm2(Xref - X_out)/N, etime
                    call print_svdvals(X_out, 'DLRA')
                end if
            else if ((if_rkad .and. .not. if_kexpm) .or. if_t_all) then
                if (io_rank()) call print_info(if_rkad, if_kexpm, Tend, dt, rk, torder)
                ! set DLRA opts
                opts = dlra_opts(mode=torder, verbose=verb, if_rank_adaptive=if_rkad, chktime=tchk)
                ! Initialize solution
                call X%initialize_LR_state(U0, S0, rk, rkmax)
                ! Run integrator
                call system_clock(count=clock_start)     ! Start Timer
                call projector_splitting_DLRA_lyapunov_integrator(X, LTI%prop, LTI%B, Tend, dt, info, & 
                       & exptA=exptA_rklib, options=opts)
                call system_clock(count=clock_stop)      ! Stop Timer
                etime = real(clock_stop-clock_start)/real(clock_rate)
                ! Reconstruct solution
                call get_state(U_out(:,:rk), X%U(:rk))
                X_out = matmul(U_out(:,:rk), matmul(X%S(:rk,:rk), transpose(U_out(:,:rk))))
                ! Print info
                if (io_rank()) then
                    write(*,'(A17,I4," TO",I1,F10.6,F8.4,E20.8,E20.8,F18.4," s")') 'OUTPUT_X   DLRA', &
                        & rk, torder, dt, Tend, norm2(X_RKlib_ref - X_out)/N, norm2(Xref - X_out)/N, etime
                    call print_svdvals(X_out, 'DLRA')
                end if
            else if ((if_rkad .and. if_kexpm) .or. if_t_all) then
                if (io_rank()) call print_info(if_rkad, if_kexpm, Tend, dt, rk, torder)
                ! set DLRA opts
                opts = dlra_opts(mode=torder, verbose=verb, if_rank_adaptive=if_rkad, chktime=tchk)
                ! Initialize solution
                call X%initialize_LR_state(U0, S0, rk, rkmax)
                ! Run integrator
                call system_clock(count=clock_start)     ! Start Timer
                call projector_splitting_DLRA_lyapunov_integrator(X, LTI%A,    LTI%B, Tend, dt, info, &
                       & exptA=exptA,       options=opts)
                call system_clock(count=clock_stop)      ! Stop Timer
                etime = real(clock_stop-clock_start)/real(clock_rate)
                ! Reconstruct solution
                call get_state(U_out(:,:rk), X%U(:rk))
                X_out = matmul(U_out(:,:rk), matmul(X%S(:rk,:rk), transpose(U_out(:,:rk))))
                ! print info
                if (io_rank()) then
                    write(*,'(A17,I4," TO",I1,F10.6,F8.4,E20.8,E20.8,F18.4," s")') 'OUTPUT_X   DLRA', &
                        & rk, torder, dt, Tend, norm2(X_RKlib_ref - X_out)/N, norm2(Xref - X_out)/N, etime
                    call print_svdvals(X_out, 'DLRA')
                end if
            end if
         end do

         if (if_save .and. io_rank()) then
            write(oname,'("example/DLRA_laplacian2D/data_X_DRLA_TO",I1,"_rk",I2.2,".npy")') torder, rk
            call save_npy(oname, X_out)
         end if

      end do
   end do

   if (if_save .and. io_rank()) then
      call save_npy("example/DLRA_laplacian2D/data_A.npy", Adata)
      call save_npy("example/DLRA_laplacian2D/data_X0.npy", X0)
      call save_npy("example/DLRA_laplacian2D/data_Q.npy", BBTdata)
      call save_npy("example/DLRA_laplacian2D/data_X_ref.npy", Xref)
      call save_npy("example/DLRA_laplacian2D/data_X_RK.npy", X_RKlib_ref)
   end if

end program demo
