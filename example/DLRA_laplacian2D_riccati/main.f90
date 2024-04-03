program demo
   use LightKrylov
   use LightKrylov_expmlib
   use LightKrylov_utils

   use LightROM_AbstractLTIsystems
   use LightROM_LyapunovSolvers
   use LightROM_LyapunovUtils
   use LightROM_RiccatiSolvers
   use LightROM_utils

   use Laplacian2D_LTI
   use Laplacian2D_LTI_Lyapunov
   use Laplacian2D_LTI_Riccati

   use stdlib_optval, only : optval 
   use stdlib_linalg, only : eye
   use stdlib_math, only : all_close, logspace
   use stdlib_io_npy, only : save_npy
   implicit none

   !----------------------------------------------------------
   !-----     LYAPUNOV EQUATION FOR LAPLACE OPERATOR     -----
   !----------------------------------------------------------

   ! DLRA
   integer, parameter :: rkmax = 14
   integer, parameter :: rk_X0 = 10
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
   type(rklib_lyapunov_mat), allocatable :: RK_prop_lyap
   type(rklib_riccati_mat), allocatable  :: RK_prop_ricc

   ! LTI system
   type(lti_system)                :: LTI
   type(state_vector), allocatable :: CT(:)
   real(kind=wp), allocatable      :: D(:,:)
   integer                         :: p
   type(LR_Q)                      :: Q

   ! Laplacian
   type(laplace_operator),   allocatable :: A

   ! LR representation
   type(state_vector), allocatable :: U(:)
   real(wp) , allocatable          :: S(:,:)
   
   !> STATE MATRIX (RKlib)
   type(state_matrix)              :: X_mat_RKlib(2)
   real(wp), allocatable           :: X_RKlib(:,:,:)
   real(wp)                        :: X_RKlib_ref(N,N)

    ! Initial condition
   real(wp)                        :: U0(N, rkmax)
   real(wp)                        :: S0(rkmax,rkmax)
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
   real(wp)  :: CTdata(N,rkmax)
   real(wp)  :: CTQcCdata(N,N)
   real(wp)  :: Rinv(rkmax,rkmax)
   real(wp)  :: Bwrk(N,rkmax)

   ! LAPACK SOLUTION LYAPUNOV
   real(wp)  :: Xref(N,N)
   ! DSYTD2
   real(wp)  :: Dm(N), work(N), wr(N), wi(N)
   real(wp)  :: E(N-1), tw(N-1)
   real(wp)  :: T(N,N), Q(N,N), Z(N,N), Vdata(N,N), Wdata(N,N), Ydata(N,N)
   real(wp)  :: scale
   integer   :: isgn
   ! SVD
   real(wp)  :: U_svd(N,N)
   real(wp)  :: S_svd(rkmax)
   real(wp)  :: V_svd(rkmax,rkmax)

   ! LAPACK SOLUTION RICATTI
   real(kind=wp)      :: Hdata(2*N,2*N)
   real(kind=wp)      :: wr2(2*N), wi2(2*N)
   real(kind=wp)      :: VR(2*N,2*N)
   integer, parameter :: Hlwork = 1040
   real(kind=wp)      :: Hwork(Hlwork)
   real(kind=wp)      :: UR(2*N,N)
   real(kind=wp)      :: UI(2*N,N)
   logical            :: flag
   real(kind=wp)      :: F(N,N)
   real(kind=wp)      :: Ginv(N,N)
   real(kind=wp)      :: Xref_ric(N,N)
   integer :: icnt

   ! timer
   integer   :: clock_rate, clock_start, clock_stop
   _
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

   p = 1
   LTI = lti_system()
   allocate(LTI%A,         source=A)
   allocate(LTI%B(1:rk_b), source=B(1:rk_b));
   allocate(LTI%CT(1:p),   source=B(1)); call mat_zero(LTI%CT)
   allocate(LTI%D(1:p,1:rk_b)); LTI%D = 0.0_wp

   ! Define initial condition
   call random_number(U0(:, 1:rk_X0))
   ! Compute SVD to get low-rank representation
   call svd(U0(:,1:rk_X0), U_svd(:,1:N), S_svd(1:rk_X0), V_svd(1:rk_X0,1:rk_X0))
   S0 = 0.0_wp
   do i = 1,rk_X0
      S0(i,i) = S_svd(i)
   end do
   U0(:,1:rk_X0) = U_svd(:,1:rk_X0)

   ! Compute the full initial condition X0 = U_in @ S0 @ U_in.T
   X0 = matmul( U0(:,1:rk_X0), matmul(S0(1:rk_X0,1:rk_X0), transpose(U0(:,1:rk_X0))))
   
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

   !------------------
   ! COMPUTE SOLUTION WITH RK FOR DIFFERENT INTEGRATION TIMES AND COMPARE TO STUART-BARTELS
   !------------------

   write(*,*) 'II.  Compute solution of the differential Lyapunov equation (Runge-Kutta):'
   write(*,*)
   ! initialize exponential propagator
   nrep = 10
   Tend = 0.1_wp
   RK_prop_lyap = rklib_lyapunov_mat(Tend)

   allocate(X_RKlib(N, N, nrep))
   call set_state(X_mat_RKlib(1:1), X0)
   write(*,'(A10,A26,A26,A20)') 'RKlib:','Tend','|| X_RK - X_ref ||_2/N', 'Elapsed time'
   write(*,*) '         ------------------------------------------------------------------------'
   do irep = 1, nrep
      call system_clock(count=clock_start)     ! Start Timer
      ! integrate
      call RK_prop_lyap%matvec(X_mat_RKlib(1), X_mat_RKlib(2))
      ! recover output
      call get_state(X_RKlib(:,:,irep), X_mat_RKlib(2:2))
      ! replace input
      call set_state(X_mat_RKlib(1:1), X_RKlib(:,:,irep))
      call system_clock(count=clock_stop)      ! Stop Timer
      write(*,'(I10,F26.4,E26.8,F18.4," s")') irep, irep*Tend, &
                     & norm2(X_RKlib(:,:,irep) - Xref)/N, &
                     & real(clock_stop-clock_start)/real(clock_rate)
   end do
!   
!   !------------------
!   ! COMPUTE DLRA FOR SHORTEST INTEGRATION TIMES FOR DIFFERENT DT AND COMPARE WITH RK SOLUTION
!   !------------------
!
!   write(*,*)
!   write(*,*) 'III. Compute approximate solution of the differential Lyapunov equation using DLRA:'
!   write(*,*)
!   write(*,'(A10,A4,A4,A10,A8,A26,A20)') 'DLRA:','  rk',' TO','dt','Tend','|| X_DLRA - X_RK ||_2/N', 'Elapsed time'
!   write(*,*) '         ------------------------------------------------------------------------'
!
!   ! Choose relevant reference case from RKlib
!   X_RKlib_ref = X_RKlib(:,:,1)
!
!   ! Choose input ranks and integration steps
!   nrk = 4; ndt = 5
!   allocate(rkv(1:nrk)); allocate(dtv(1:ndt)); 
!   rkv = (/ 2, 6, 10, 14 /)
!   dtv = logspace(-5.0_wp, -1.0_wp, ndt)
!
!   do torder = 1, 2
!      do i = 1, nrk
!         rk = rkv(i)
!
!         allocate(U(1:rk)); call mat_zero(U)
!         allocate(S(1:rk,1:rk)); S = 0.0_wp
!         write(*,'(A10,I1)') ' torder = ', torder
!
!         do j = ndt, 1, -1
!            dt = dtv(j)
!            if (verb) write(*,*) '    dt = ', dt, 'Tend = ', Tend
!
!            ! Reset input
!            call set_state(U(1:rk), U0(:,1:rk))
!            S(1:rk,1:rk) = S0(1:rk,1:rk)
!
!            ! run step
!            call system_clock(count=clock_start)     ! Start Timer
!            call numerical_low_rank_splitting_integrator(U(1:rk), S(1:rk,1:rk), &
!                                                       & LTI, Tend, dt, torder, info)
!            call system_clock(count=clock_stop)      ! Stop Timer
!
!            ! Reconstruct solution
!            call get_state(U_out(:,1:rk), U)
!            X_out = matmul(U_out(:,1:rk), matmul(S(1:rk,1:rk), transpose(U_out(:,1:rk))))
!      
!            write(*,'(A10,I4," TO",I1,F10.6,F8.4,E26.8,F18.4," s")') 'OUTPUT', &
!                              & rk, torder, dt, Tend, &
!                              & norm2(X_RKlib_ref - X_out)/N, &
!                              & real(clock_stop-clock_start)/real(clock_rate)
!         end do
!
!         if (save) then
!            write(oname,'("example/DLRA_laplacian2D/data_X_DRLA_TO",I1,"_rk",I2.2,".npy")'), torder, rk
!            call save_npy(oname, X_out)
!         end if
!
!         deallocate(U);
!         deallocate(S);
!
!      end do
!   end do
!
!   if (save) then
!      call save_npy("example/DLRA_laplacian2D/data_A.npy", Adata)
!      call save_npy("example/DLRA_laplacian2D/data_X0.npy", X0)
!      call save_npy("example/DLRA_laplacian2D/data_Q.npy", BBTdata)
!      call save_npy("example/DLRA_laplacian2D/data_X_ref.npy", Xref)
!      call save_npy("example/DLRA_laplacian2D/data_X_RK.npy", X_RKlib_ref)
!   end if

   write(*,*)
   write(*,*) 'RICCATI EQUATION FOR THE 2D LAPLACE OPERATOR:'
   write(*,*)
   write(*,*) '    Algebraic Lyapunov equation:'
   write(*,*) '                0 = A @ X + X @ A.T + Q - X @ B @ R^{-1} @ B.T @ X'
   write(*,*)
   write(*,*) '    Differential Lyapunov equation:'
   write(*,*) '          \dot{X} = A @ X + X @ A.T + Q - X @ B @ R^{-1} @ B.T @ X'
   write(*,*)
   write(*,'(A16,I4,"x",I4)') '  Problem size: ', N, N
   write(*,*)
   write(*,*) '  Initial condition: rank(X0) =', rk_X0
   write(*,*) '  Inhomogeneity:     rank(Q)  =', rk_C
   write(*,*) '  Inhomogeneity:     rank(B)  =', rk_B
   write(*,*)
   write(*,*) '---------------------------------------------'
   write(*,*)

   ! Define C
   call init_rand(CT, ifnorm = .false.)
   call get_state(CTdata(:,1:rk_c), CT)
   Qcdata = eye(rk_c)
   CTQcCdata = matmul(CTdata(:,1:rk_c), matmul( Qcdata, transpose(CTdata(:,1:rk_c))))
   CTQcC(1:N**2) = reshape(CTQcCdata, shape(CTQcC))
   ! Update LTI system
   call mat_copy(LTI%CT(1:rk_c), CT)

   ! Define Rinv & compule BRinvBT
   Rinv = 1e3*eye(rkmax)
   Bwrk = matmul(Bdata, Rinv)
   call get_state(Bdata(:,1:rk_b), B)
   BRinvBTdata = matmul( Bwrk, transpose(Bdata) )

   ! construct Hmat
   Hdata = 0.0_wp
   Hdata(  1:N,    1:N  ) =  Adata
   Hdata(N+1:2*N,N+1:2*N) = -transpose(Adata)
   Hdata(  1:N,  N+1:2*N) = -BRinvBTdata
   Hdata(N+1:2*N,  1:N  ) = -CTQcCdata

   call dgeev('N', 'V', 2*N, Hdata, 2*N, wr2, wi2, VR, 2*N, VR, 2*N, Hwork, Hlwork, info)

   ! extract stable eigenspace
   UR = 0.0_wp
   UI = 0.0_wp
   icnt = 0
   flag = .true.
   do i = 1, 2*N
      if ( wr2(i) .lt. 0.0 ) then
         icnt = icnt + 1
         if ( wi2(i) .eq. 0.0 ) then ! real
            UR(:,icnt) = VR(:,i)
         else                       ! complex
            if (flag) then
               UR(:,icnt)   =  VR(:,i)
               UI(:,icnt)   =  VR(:,i+1)
               UR(:,icnt+1) =  VR(:,i)
               UI(:,icnt+1) = -VR(:,i+1)
               flag = .false.
            else
               flag = .true.
            end if
         end if
      end if
   end do
   ! construct solution
   F    = matmul(UR(n+1:2*N,:), transpose(UR(1:N,:))) + matmul(UI(n+1:2*n,:), transpose(UI(1:N,:)))
   Ginv = matmul(UR(  1:N,  :), transpose(UR(1:N,:))) + matmul(UI(  1:N,  :), transpose(UI(1:N,:)))
   call inv(Ginv)
   Xref_ric = matmul(F, Ginv)

   call print_mat(N,N,Xref_ric)

   X0 = matmul(Adata, Xref_ric) + matmul(Xref_ric, transpose(Adata)) + CTQcCdata - matmul(Xref_ric, matmul(BRinvBTdata, Xref_ric))

   call print_mat(N,N,X0)

   ! Define initial condition
   U0 = 0.0_wp
   call random_number(U0(:, 1:rk_X0))
   ! Compute SVD to get low-rank representation
   call svd(U0(:,1:rk_X0), U_svd(:,1:N), S_svd(1:rk_X0), V_svd(1:rk_X0,1:rk_X0))
   S0 = 0.0_wp
   do i = 1,rk_X0
      S0(i,i) = S_svd(i)
   end do
   U0(:,1:rk_X0) = U_svd(:,1:rk_X0)

   ! Compute the full initial condition X0 = U_in @ S0 @ U_in.T
   X0 = matmul( U0(:,1:rk_X0), matmul(S0(1:rk_X0,1:rk_X0), transpose(U0(:,1:rk_X0))))

   ! initialize exponential propagator
   nrep = 10
   Tend = 2.0_wp
   RK_prop_ricc = rklib_riccati_mat(Tend)

   allocate(X_RKlib(N, N, nrep))
   call set_state_mat(Xmat_RKlib(1:1), X0)
   write(*,'(A10,A26,A26,A20)') 'RKlib:','Tend','|| dX/dt ||_2/N', 'Elapsed time'
   write(*,*) '         ------------------------------------------------------------------------'
   do irep = 1, nrep
      call system_clock(count=clock_start)     ! Start Timer
      ! integrate
      call RK_prop_ricc%matvec(Xmat_RKlib(1), Xmat_RKlib(2))
      ! recover output
      call get_state_mat(X_RKlib(:,:,irep), Xmat_RKlib(2:2))
      ! replace input
      call set_state_mat(Xmat_RKlib(1:1), X_RKlib(:,:,irep))
      call system_clock(count=clock_stop)      ! Stop Timer
      write(*,'(I10,F26.4,E26.8,F18.4," s")') irep, irep*Tend, &
                     & norm2(CARE(X_RKlib(:,:,irep), Adata, Qdata, BRinvBTdata))/N, &
                     & real(clock_stop-clock_start)/real(clock_rate)
   end do


contains

   !--------------------------------------------------------------------
   !-----     UTILITIES FOR STATE_VECTOR AND STATE MATRIX TYPES    -----
   !--------------------------------------------------------------------

   subroutine get_state(mat_out, state_in)
      !! Utility function to transfer data from a state vector to a real array
      real(kind=wp),          intent(out) :: mat_out(:,:)
      class(abstract_vector), intent(in)  :: state_in(:)
      ! internal variables
      integer :: k, kdim
      mat_out = 0.0_wp
      select type (state_in)
      type is (state_vector)
         kdim = size(state_in)
         call assert_shape(mat_out, (/ N, kdim /), 'get_state -> state_vector', 'mat_out')
         do k = 1, kdim
            mat_out(:,k) = state_in(k)%state
         end do
      type is (state_matrix)
         call assert_shape(mat_out, (/ N, N /), 'get_state -> state_matrix', 'mat_out')
         mat_out = reshape(state_in(1)%state, (/ N, N /))
      end select
      return
   end subroutine get_state

   subroutine set_state(state_out, mat_in)
      !! Utility function to transfer data from a real array to a state vector
      class(abstract_vector), intent(out) :: state_out(:)
      real(kind=wp),          intent(in)  :: mat_in(:,:)
      ! internal variables
      integer       :: k, kdim
      select type (state_out)
      type is (state_vector)
         kdim = size(state_out)
         call assert_shape(mat_in, (/ N, kdim /), 'set_state -> state_vector', 'mat_in')
         call mat_zero(state_out)
         do k = 1, kdim
            state_out(k)%state = mat_in(:,k)
         end do
      type is (state_matrix)
         call assert_shape(mat_in, (/ N, N /), 'set_state -> state_matrix', 'mat_in')
         call mat_zero(state_out)
         state_out(1)%state = reshape(mat_in, shape(state_out(1)%state))
      end select
      return
   end subroutine set_state

   subroutine init_rand(state, ifnorm)
      !! Utility function to initialize a state vector with random data
      class(abstract_vector), intent(inout)  :: state(:)
      logical, optional,      intent(in)     :: ifnorm
      ! internal variables
      integer :: k, kdim
      logical :: normalize
      normalize = optval(ifnorm,.true.)
      select type (state)
      type is (state_vector)
         kdim = size(state)
         do k = 1, kdim
            call state(k)%rand(ifnorm = normalize)
         end do
      type is (state_matrix)
         kdim = size(state)
         do k = 1, kdim
            call state(k)%rand(ifnorm = normalize)
         end do
      end select
      return
   end subroutine init_rand

   subroutine build_operator(A)
      !! Build the two-dimensional Laplace operator explicitly
      real(kind=wp), intent(out) :: A(N,N)

      A = -4.0_wp/dx2*eye(N)
      do i = 1, nx
         do j = 1, nx - 1
            k = (i-1)*nx + j
            A(k + 1, k) = 1.0_wp/dx2
            A(k, k + 1) = 1.0_wp/dx2
         end do 
      end do
      do i = 1, N-nx
         A(i, i + nx) = 1.0_wp/dx2
         A(i + nx, i) = 1.0_wp/dx2
      end do
      return
   end subroutine build_operator

   subroutine reconstruct_TQ(T, Q, A, D, E, tw)
      !! Reconstruct tridiagonal matrix T and orthogonal projector Q from dsytd2 output (A, D, E)
      real(kind=wp), intent(out) :: T(N,N)
      real(kind=wp), intent(out) :: Q(N,N)
      real(kind=wp), intent(in)  :: A(N,N)
      real(kind=wp), intent(in)  :: D(N)
      real(kind=wp), intent(in)  :: E(N-1)
      real(kind=wp), intent(in)  :: tw(N-1)

      ! internal variables
      real(wp)  :: Hi(N,N)
      real(wp)  :: vec(N,1)

      ! Build orthogonal Q = H(1) @  H(2) @ ... @ H(n-1)
      Q = eye(N)
      do i = 1, N - 1
         vec          = 0.0_wp
         vec(i+1,1)   = 1.0_wp
         vec(i+2:N,1) = A(i+2:N,i)
         Hi           = eye(N) - tw(i) * matmul( vec, transpose(vec) )
         Q            = matmul( Q, Hi )
      end do

      ! Build tridiagonal T
      T = 0.0_wp
      do i = 1, N
         T(i,i) = D(i)
      end do
      do i = 1, N - 1
         T(i,i+1) = E(i)
         T(i+1,i) = E(i)
      end do

   end subroutine reconstruct_TQ

end program demo