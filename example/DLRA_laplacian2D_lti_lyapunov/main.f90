program demo
   ! Standard Library
   use stdlib_optval, only : optval 
   use stdlib_linalg, only : eye
   use stdlib_math, only : all_close, logspace
   use stdlib_io_npy, only : save_npy
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
   integer, parameter :: rkmax = 14
   integer, parameter :: rk_X0 = 14
   logical, parameter :: verb  = .false.
   logical, parameter :: save  = .false.
   character*128      :: oname
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
   integer                         :: i, j, k, irep, nrep

   ! PROBLEM DEFINITION
   real(wp)  :: Adata(N,N)
   real(wp)  :: Bdata(N,rkmax)
   real(wp)  :: BBTdata(N,N)
   
   ! LAPACK SOLUTION
   real(wp)  :: Xref(N,N)
   ! DSYTD2
   real(wp)  :: Dm(N), work(N), wr(N), wi(N)
   real(wp)  :: E(N-1), tw(N-1)
   real(wp)  :: T(N,N), Q(N,N), Z(N,N), Vdata(N,N), Wdata(N,N), Ydata(N,N)
   real(wp)  :: scale
   integer   :: isgn

   ! timer
   integer   :: clock_rate, clock_start, clock_stop

   call logger%configure(level=error_level); write(*,*) 'Logging set to error_level.'

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
   call get_state(Bdata(:,1:rk_b), B)
   BBTdata = -matmul(Bdata(:,1:rk_b), transpose(Bdata(:,1:rk_b)))
   BBT(1:N**2) = -reshape(BBTdata, shape(BBT))

   ! Define LTI system
   LTI = lti_system()
   call LTI%initialize_lti_system(A, B, B)
   call zero_basis(LTI%CT)

   ! Define initial condition of the form X0 + U0 @ S0 @ U0.T SPD 
   if (verb) write(*,*) '    Define initial condition'
   call generate_random_initial_condition(U0, S0, rk_X0)
   call get_state(U_out, U0)
   
   ! Compute the full initial condition X0 = U_in @ S0 @ U_in.T
   X0 = matmul( U_out(:,1:rk_X0), matmul(S0(1:rk_X0,1:rk_X0), transpose(U_out(:,1:rk_X0))))
   
   !------------------
   ! COMPUTE EXACT SOLUTION OF THE LYAPUNOV EQUATION WITH LAPACK
   !------------------

   write(*,*) 'I.   Exact solution of the algebraic Lyapunov equation (LAPACK):'
   call system_clock(count=clock_start)     ! Start Timer

   ! Explicit 2D laplacian
   call build_operator(Adata)

   ! Transform operator to tridiagonal form
   call dsytd2('L', N, Adata, N, Dm, E, tw, info)
   
   ! Reconstruct T and Q
   call reconstruct_TQ(T, Q, Adata(1:N,1:N), Dm, E, tw)
   
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
   write(*,*)

   ! sanity check
   X0 = CALE(Xref, Adata, BBT)
   write(*,*) '    Direct problem:', norm2(X0)/N

   !------------------
   ! COMPUTE SOLUTION WITH RK FOR DIFFERENT INTEGRATION TIMES AND COMPARE TO STUART-BARTELS
   !------------------

   write(*,*) 'II.  Compute solution of the differential Lyapunov equation (Runge-Kutta):'
   write(*,*)
   ! initialize exponential propagator
   nrep = 10
   Tend = 0.1_wp
   RK_propagator = rklib_lyapunov_mat(Tend)

   allocate(X_RKlib(N, N, nrep))
   call set_state(X_mat_RKlib(1:1), X0)
   write(*,'(A10,A26,A26,A20)') 'RKlib:','Tend','|| X_RK - X_ref ||_2/N', 'Elapsed time'
   write(*,*) '         ------------------------------------------------------------------------'
   do irep = 1, nrep
      call system_clock(count=clock_start)     ! Start Timer
      ! integrate
      call RK_propagator%matvec(X_mat_RKlib(1), X_mat_RKlib(2))
      ! recover output
      call get_state(X_RKlib(:,:,irep), X_mat_RKlib(2:2))
      ! replace input
      call set_state(X_mat_RKlib(1:1), X_RKlib(:,:,irep))
      call system_clock(count=clock_stop)      ! Stop Timer
      write(*,'(I10,F26.4,E26.8,F18.4," s")') irep, irep*Tend, &
                     & norm2(X_RKlib(:,:,irep) - Xref)/N, &
                     & real(clock_stop-clock_start)/real(clock_rate)
   end do
   
   !------------------
   ! COMPUTE DLRA FOR SHORTEST INTEGRATION TIMES FOR DIFFERENT DT AND COMPARE WITH RK SOLUTION
   !------------------

   write(*,*)
   write(*,*) 'III. Compute approximate solution of the differential Lyapunov equation using DLRA:'
   write(*,*)
   write(*,'(A10,A4,A4,A10,A8,A26,A20)') 'DLRA:','  rk',' TO','dt','Tend','|| X_DLRA - X_RK ||_2/N', 'Elapsed time'
   write(*,*) '         ------------------------------------------------------------------------'

   ! Choose relevant reference case from RKlib
   X_RKlib_ref = X_RKlib(:,:,1)

   ! Choose input ranks and integration steps
   nrk = 4; ndt = 5
   allocate(rkv(1:nrk)); allocate(dtv(1:ndt)); 
   rkv = (/ 2, 6, 10, 14 /)
   dtv = logspace(-5.0_wp, -1.0_wp, ndt)

   X = LR_state()

   do torder = 1, 2
      do i = 1, nrk
         rk = rkv(i)

         write(*,'(A10,I1)') ' torder = ', torder

         do j = ndt, 1, -1
            dt = dtv(j)
            if (verb) write(*,*) '    dt = ', dt, 'Tend = ', Tend

            ! Reset input
            call X%initialize_LR_state(U0, S0, rk)

            ! run step
            call system_clock(count=clock_start)     ! Start Timer
            call numerical_low_rank_splitting_lyapunov_integrator(X, LTI%A, LTI%B, Tend, dt, torder, info)
            call system_clock(count=clock_stop)      ! Stop Timer

            ! Reconstruct solution
            call get_state(U_out(:,1:rk), X%U)
            X_out = matmul(U_out(:,1:rk), matmul(X%S, transpose(U_out(:,1:rk))))
      
            write(*,'(A10,I4," TO",I1,F10.6,F8.4,E26.8,F18.4," s")') 'OUTPUT', &
                              & rk, torder, dt, Tend, &
                              & norm2(X_RKlib_ref - X_out)/N, &
                              & real(clock_stop-clock_start)/real(clock_rate)

            deallocate(X%U)
            deallocate(X%S)
         end do

         if (save) then
            write(oname,'("example/DLRA_laplacian2D/data_X_DRLA_TO",I1,"_rk",I2.2,".npy")') torder, rk
            call save_npy(oname, X_out)
         end if

      end do
   end do

   if (save) then
      call save_npy("example/DLRA_laplacian2D/data_A.npy", Adata)
      call save_npy("example/DLRA_laplacian2D/data_X0.npy", X0)
      call save_npy("example/DLRA_laplacian2D/data_Q.npy", BBTdata)
      call save_npy("example/DLRA_laplacian2D/data_X_ref.npy", Xref)
      call save_npy("example/DLRA_laplacian2D/data_X_RK.npy", X_RKlib_ref)
   end if

end program demo