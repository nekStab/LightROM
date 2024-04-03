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
   type(rklib_riccati_mat), allocatable  :: RK_prop_ricc

   ! LTI system
   type(lti_system)                :: LTI
   real(kind=wp), allocatable      :: D(:,:)
   integer                         :: p

   ! Laplacian
   type(laplace_operator),   allocatable :: A

   ! LR representation
   type(state_vector), allocatable :: U(:)
   real(wp) , allocatable          :: S(:,:)
   
   !> STATE MATRIX (RKlib)
   type(state_matrix)              :: Xmat_RKlib(2)
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
   real(wp)  :: CTdata(N,rkmax)
   real(wp)  :: Qc(rkmax,rkmax)
   real(wp)  :: Rinv(rkmax,rkmax)
   real(wp)  :: Bwrk(N,rkmax)

   ! SVD
   real(wp)  :: U_svd(N,N)
   real(wp)  :: S_svd(rkmax)
   real(wp)  :: V_svd(rkmax,rkmax)

   ! LAPACK SOLUTION RICATTI
   real(kind=wp)      :: Hdata(2*N,2*N)
   real(kind=wp)      :: wr(2*N), wi(2*N)
   real(kind=wp)      :: VR(2*N,2*N)
   integer, parameter :: lwork = 1040
   real(kind=wp)      :: work(lwork)
   real(kind=wp)      :: UR(2*N,N)
   real(kind=wp)      :: UI(2*N,N)
   logical            :: flag
   real(kind=wp)      :: F(N,N)
   real(kind=wp)      :: Ginv(N,N)
   real(kind=wp)      :: Xref(N,N)
   integer :: icnt

   real(kind=wp)      :: test(N**2)
   real(kind=wp)      :: test2(N**2)

   ! timer
   integer   :: clock_rate, clock_start, clock_stop
   
   call system_clock(count_rate=clock_rate)

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

   ! Define C, Qc & compute CTQcC
   Qc = eye(rkmax)

   call init_rand(CT, ifnorm = .true.)
   call get_state(CTdata(:,1:rk_c), CT)
   CTQcCdata = matmul(CTdata(:,1:rk_c), matmul( Qc(1:rk_c,1:rk_c), transpose(CTdata(:,1:rk_c))))
   CTQcC(1:N**2) = reshape(CTQcCdata, shape(CTQcC))
   
   ! Define B, Rinv & compule BRinvBT
   Rinv = 1e1*eye(rkmax)

   call init_rand(B, ifnorm = .true.)
   Bdata = 0.0_wp
   call get_state(Bdata(:,1:rk_b), B)
   Bwrk = 0.0_wp
   Bwrk(:,1:rk_b) = matmul(Bdata(:,1:rk_b), Rinv(1:rk_b, 1:rk_b))
   BRinvBTdata = matmul( Bwrk(:,1:rk_b), transpose(Bdata(:,1:rk_b)) )

   ! Define LTI system
   LTI = lti_system()
   allocate(LTI%A,          source=A)
   allocate(LTI%B(1:rk_b),  source=B(1:rk_b))
   allocate(LTI%CT(1:rk_c), source=CT(1:rk_c))
   allocate(LTI%D(1:p,1:rk_b)); LTI%D = 0.0_wp

   !------------------
   ! COMPUTE EXACT SOLUTION OF THE RICCATI EQUATION WITH LAPACK
   !------------------

   write(*,*) 'I.   Exact solution of the algebraic Riccati equation (LAPACK):'
   call system_clock(count=clock_start)     ! Start Timer

   ! Explicit 2D laplacian
   call build_operator(Adata)

   ! construct Hmat
   Hdata = 0.0_wp
   Hdata(  1:N,    1:N  ) =  Adata
   Hdata(N+1:2*N,N+1:2*N) = -transpose(Adata)
   Hdata(  1:N,  N+1:2*N) = -BRinvBTdata
   Hdata(N+1:2*N,  1:N  ) = -CTQcCdata

   call dgeev('N', 'V', 2*N, Hdata, 2*N, wr, wi, VR, 2*N, VR, 2*N, work, lwork, info)

   ! extract stable eigenspace
   UR = 0.0_wp
   UI = 0.0_wp
   icnt = 0
   flag = .true.
   do i = 1, 2*N
      if ( wr(i) .lt. 0.0 ) then
         icnt = icnt + 1
         if ( wi(i) .eq. 0.0 ) then ! real
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
   Xref = matmul(F, Ginv)

   ! sanity check
   !call print_mat(N,N,Xref)
   !X0 = CARE(Xref, Adata, CTQcCdata, BRinvBTdata)

   call system_clock(count=clock_stop)      ! Stop Timer
   write(*,'(A40,F10.4," s")') '--> X_ref.    Elapsed time:', real(clock_stop-clock_start)/real(clock_rate)
   write(*,*)

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
   !nrep = 10
   !Tend = 0.1_wp !2.0_wp
   !RK_prop_ricc = rklib_riccati_mat(Tend)
!
   !allocate(X_RKlib(N, N, nrep))
   !call set_state(Xmat_RKlib(1:1), X0)
   !write(*,'(A10,A26,A26,A20)') 'RKlib:','Tend','|| dX/dt ||_2/N', 'Elapsed time'
   !write(*,*) '         ------------------------------------------------------------------------'
   !do irep = 1, nrep
   !   call system_clock(count=clock_start)     ! Start Timer
   !   ! integrate
   !   call RK_prop_ricc%matvec(Xmat_RKlib(1), Xmat_RKlib(2))
   !   ! recover output
   !   call get_state(X_RKlib(:,:,irep), Xmat_RKlib(2:2))
   !   call print_mat(N**2,1,Xmat_RKlib(2)%state)
   !   ! replace input
   !   call set_state(Xmat_RKlib(1:1), X_RKlib(:,:,irep))
   !   call system_clock(count=clock_stop)      ! Stop Timer
   !   write(*,'(I10,F26.4,E26.8,F18.4," s")') irep, irep*Tend, &
   !                  & norm2(CARE(X_RKlib(:,:,irep), Adata, CTQcCdata, BRinvBTdata))/N, &
   !                  & real(clock_stop-clock_start)/real(clock_rate)
   !end do



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

end program demo