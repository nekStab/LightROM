program demo
   use LightKrylov
   use Laplacian
   use Laplacian_Lyapunov
   use LightROM_LyapunovSolvers
   use LightROM_LyapunovUtils
   use LightROM_utils
   use LightROM_expmlib
   use stdlib_linalg, only : eye
   use stdlib_math, only : all_close
   !use stdlib_io_npy, only : save_npy
   implicit none

   !------------------------------------------------
   !-----     LINEAR OPERATOR INVESTIGATED     -----
   !------------------------------------------------

   ! Laplacian in Lyapunov equation
   !> Exponential propagator.
   type(lyapunov_operator),  allocatable :: A_lyap
   type(rklib_lyapunov_mat), allocatable :: A_prop_mat
   !> Krylov subspace.
   type(state_matrix), allocatable :: X_mat(:)
   type(state_matrix), allocatable :: wrk_mat(:)

   ! Regular Laplacian
   !> Exponential propagator.
   type(laplace_operator),           allocatable :: A
   type(rklib_exptA_laplacian_vec),  allocatable :: A_prop_vec
   type(krylov_exptA_laplacian_vec), allocatable :: Ak_prop_vec
   !> Krylov subspace.
   type(state_vector), allocatable :: X_vec(:)
   type(state_vector), allocatable :: wrk_vec(:)

   !> Sampling time.
   real(kind=wp), parameter :: tau = 0.1_wp !0.1_wp
   real(kind=wp), parameter :: tau_long = 1.0_wp !1.0_wp

   !> Number of eigenvalues we wish to converge.
   integer, parameter :: nev = 10
   !> Krylov subspace dimension.
   integer, parameter :: kdim = 2*nx
   !> Coordinates of the eigenvectors in the Krylov basis.
   complex(kind=wp) :: v(kdim, kdim)
   !> Eigenvalues.
   complex(kind=wp) :: lambda(kdim)
   !> Residual.
   real(kind=wp)    :: residuals(kdim)
   !> Information flag.
   integer          :: info
   !> Misc
   integer       :: i, j, k
   real(kind=wp) :: alpha
   class(abstract_vector), allocatable :: wrk
   complex(kind=wp)                    :: eigenvectors(nx, nev)
   real(wp)      :: tmp(nx,nx)
   !> DLRA
   integer, parameter :: rk = 5
   integer :: torder
   integer :: ipiv(nx)
   real(wp) :: Tend
   real(wp) :: S(rk,rk)
   real(wp) :: dt          !> timestep

   !> lapack
   real(wp)  :: Amat(nx,nx)
   real(wp)  :: Xmat(nx,nx)
   real(wp)  :: Emat(nx,nx)
   real(wp)  :: exptA(nx,nx)
   real(wp)  :: chkmat(nx,nx)
   real(wp)  :: work(nx)
   real(wp)  :: wr(nx)
   real(wp)  :: wi(nx)
   real(wp)  :: Bmat(nx,rk)
   real(wp)  :: Qmat(nx,nx)
   real(wp)  :: Wmat(nx,nx)
   real(wp)  :: T(nx,nx)
   real(wp)  :: Z(nx,nx)
   real(wp)  :: scale
   integer   :: isgn
   !> testing DLRA
   real(wp) :: zz,oo
   ! Mstep
   real(wp)  :: tol
   logical   :: verb
   real(wp)  :: R(rk,rk)
   real(wp)  :: Swrk(rk,rk)
   real(wp)  :: S0(rk,rk)
   real(wp)  :: X0(nx,nx)
   integer   :: nsteps
   ! G step
   !     K, S, L steps
   type(state_vector) :: K0(rk)
   type(state_vector) :: Kdot(rk)
   real(wp) :: Sh(rk,rk)
   real(wp) :: Sdot(rk,rk)
   type(state_vector) :: L0T(rk)
   type(state_vector) :: LdotT(rk)
   
   type(state_vector) :: U1(rk)
   type(state_vector) :: Uwrk(rk)
   real(wp)           :: mval(rk)
   ! RKlib
   type(state_matrix) :: X_RKlib_mat(2)
   real(wp)           :: Xtmp(nx**2,1)
   real(wp)           :: Xmat_rklib_out(nx, nx)

   real(wp)  :: Xmat_lyap(nx**2,3)
   real(wp)  :: Xmat_DLRA_0(nx,rk)
   real(wp)  :: Xmat_DLRA_in(nx,rk)
   real(wp)  :: Xmat_DLRA_out(nx,rk)
   type(state_vector), allocatable :: B(:)
   type(state_vector), allocatable :: U(:)
   real(wp)  :: US(nx,rk)
   real(wp)  :: USUT(nx,nx)

   !------------------
   ! DLRA
   !------------------

   allocate(B(1:rk)); call mat_zero(B)
   !do i = 1, rk
   !   call random_number(B(i)%state)
   !end do
   B(1)%state = 1.0_wp

   !------------------
   ! COMPUTE EXACT SOLUTION TO LYAPUNOV EQUATION WITH LAPACK
   !------------------

   ! Explicit laplacian
   Amat = -2.0/dx**2*eye(nx)
   do i = 1, nx-1
      Amat(i+1,i) = 1.0/dx**2
      Amat(i,i+1) = 1.0/dx**2
   end do
   ! laplace operator already in tridiagonal form
   ! compute real Schur form A = Z @ T @ Z.T
   T = Amat
   call dhseqr('S', 'I', nx, 1, nx, T, nx, wr, wi, Z, nx, work, nx, info )
   Emat = 0.0_wp
   Emat = matmul(Z, matmul(T, transpose(Z)))

   ! Explicit RHS
   Bmat = 0.0_wp
   do i = 1, rk
      Bmat(:,i) = B(i)%state
   end do
   Qmat = matmul(Bmat,transpose(Bmat))
   !write(*,*) 'Bmat'
   !call print_mat(nx, rk, Bmat)
   ! Change to Z basis
   Wmat = -matmul(transpose(Z), matmul(Qmat, Z))
   
   ! Compute solution of Lyapunov equation For Schur decomposition
   isgn = 1
   scale = 1.0_wp
   call dtrsyl('N', 'T', isgn, nx, nx, T, nx, T, nx, Wmat, nx, scale, info)
   ! Change back to original basis
   Xmat = matmul(Z, matmul(Wmat, transpose(Z)))
   !> sanity check
   write(*,*) 'LAPACK DTRSYL      A @ X + X @ A.T + B @ B.T == 0 : ',&
         & all_close(matmul(Amat,Xmat) + matmul(Xmat,transpose(Amat)), -Qmat, rtol, atol)

   !------------------
   ! COMPUTE THE MATRIX EXPONENTIAL OF THE OPERATOR
   !------------------

   !> Matrix exponential: exp(dt*A) = Z @ diag(exp(dt*T_ii)) @ Z.T
   dt = 0.1_wp
   Emat = 0.0_wp
   do i = 1, nx
      Emat(i,i) = exp(dt * T(i,i))
   end do
   exptA = matmul(Z, matmul(Emat, transpose(Z)))
   !> sanity check
   call expm(Emat, dt*Amat)
   write(*,*) 'MATRIX EXPONENTIAL                     expm(dt*A) : ', &
         & all_close(Emat, exptA, rtol, atol)

   !------------------
   ! COMPUTE SINGLE M STEP USING KEXPM AND COMPARE TO MATRIX EXPONENTIAL
   !------------------

   !> Random initial condition of rank 4
   Xmat_DLRA_in = 0.0_wp
   allocate(U(1:rk));
   do i = 1,rk
      call random_number(U(i)%state)
      Xmat_DLRA_0(:,i) = U(i)%state
   enddo
   !call mat_zero(U)
   !zz = 0.0_wp
   !oo = 1.0_wp
   !!U(1)%state = (/ oo, zz, zz, zz, zz, zz, zz, zz  /)
   !U(1)%state = (/ zz, oo, zz, zz, zz, zz, zz, zz  /)
   !U(2)%state = (/ zz, zz, oo, zz, zz, zz, zz, zz  /)
   !U(3)%state = (/ zz, zz, zz, oo, zz, zz, zz, zz  /)
   !U(4)%state = (/ zz, zz, zz, zz, oo, zz, zz, zz  /)
   !U(5)%state = (/ zz, zz, zz, zz, zz, oo, zz, zz  /)
   !do i = 1,rk
   !   Xmat_DLRA_0(:,i) = U(i)%state
   !enddo
   !> Orthonormalize (in-place)
   S = 0.0_wp
   call qr_factorization(U,S,info)
   do i = 1,rk
      Xmat_DLRA_in(:,i) = U(i)%state
   end do
   S0 = S
   X0 = matmul( Xmat_DLRA_in, matmul(S0, transpose(Xmat_DLRA_in) ) )

   ! do M step
   tol = 1e-10
   verb = .false.
   do i = 1, rk
      call kexpm(Uwrk(1), A, U(i), dt, tol, info, verbosity = verb)
      call U(i)%axpby(0.0_wp, Uwrk(1), 1.0_wp) ! overwrite old solution
   enddo
   ! Reorthonormalize
   R = 0.0_wp; Swrk = 0.0_wp
   call qr_factorization(U, R, info)
   !> Update low-rank coefficient matrix
   Swrk = matmul(S, transpose(R))
   S    = matmul(R, Swrk)

   !> check result
   do i = 1, rk
      Xmat_DLRA_out(:,i) = U(i)%state
   end do
   USUT = matmul( Xmat_DLRA_out, matmul(S, transpose(Xmat_DLRA_out) ) )
   chkmat = matmul( exptA, matmul( X0, transpose(exptA) ) )

   write(*,*) 'SINGLE MSTEP                                      : ', &
            & all_close(USUT, chkmat, rtol, atol)
   !call print_mat(nx, nx, USUT)

   !------------------
   ! COMPUTE SEVERAL M STEPS USING KEXPM AND COMPARE TO MATRIX EXPONENTIAL
   !------------------

   ! do several M steps
   nsteps = 200
   dt = 0.01_wp
   ! Reset input
   do i = 1,rk
      U(i)%state = Xmat_DLRA_in(:,i)
   end do
   S = S0
   ! do M step
   tol = 1e-10
   verb = .false.
   do j = 1, nsteps
      !write(*,*) 'M step', j
      do i = 1, rk
         !write(*,*) '    kexpm',i
         call kexpm(Uwrk(1), A, U(i), dt, tol, info, verbosity = verb)
         call U(i)%axpby(0.0_wp, Uwrk(1), 1.0_wp) ! overwrite old solution
      enddo
      ! Reorthonormalize
      R = 0.0_wp; Swrk = 0.0_wp
      call qr_factorization(U, R, info)
      !> Update low-rank coefficient matrix
      Swrk = matmul(S, transpose(R))
      S    = matmul(R, Swrk)
   end do
   ! Recompute explicit matrix exponential
   Tend = nsteps*dt
   Emat = 0.0_wp
   do i = 1, nx
      Emat(i,i) = exp(Tend * T(i,i))
   end do
   exptA = matmul(Z, matmul(Emat, transpose(Z)))

   !> check result
   do i = 1, rk
      Xmat_DLRA_out(:,i) = U(i)%state
   end do
   USUT = matmul( Xmat_DLRA_out, matmul(S, transpose(Xmat_DLRA_out) ) )
   chkmat = matmul( exptA, matmul( X0, transpose(exptA) ) )

   write(*,*) 'MULTIPLE MSTEP                                    : ', &
            & all_close(USUT, chkmat, rtol, atol)
   !call print_mat(nx, nx, USUT)

   !------------------
   ! COMPUTE SINGLE M STEP + G STEP USING KEXPM AND COMPARE TO RKLIB
   !------------------
   dt = 0.1_wp

   ! Reset input
   do i = 1,rk
      U(i)%state = Xmat_DLRA_in(:,i)
   end do
   S = S0

   ! M step
   tol = 1e-10
   verb = .false.
   do i = 1, rk
      call kexpm(Uwrk(1), A, U(i), dt, tol, info, verbosity = verb)
      call U(i)%axpby(0.0_wp, Uwrk(1), 1.0_wp) ! overwrite old solution
   enddo
   ! Reorthonormalize
   R = 0.0_wp; Swrk = 0.0_wp
   call qr_factorization(U, R, info)
   !> Update low-rank coefficient matrix
   Swrk = matmul(S, transpose(R))
   S    = matmul(R, Swrk)
   
   ! G step
   call mat_zero(U1); call mat_zero(Uwrk);


   !     K step
   call mat_mult(U1, U, S)               ! K0
   call apply_outerproduct(Uwrk, B, U)   ! Kdot
   !> Construct solution U1
   call mat_axpby(U1, 1.0_wp, Uwrk, dt)  ! K0 + tau*Kdot
   !> Orthonormalize in-place
   Swrk = 0.0_wp
   call qr_factorization(U1, Swrk, info)
   !write(*,*) 'S'
   !call print_mat(rk, rk, Swrk)
   S = Swrk


   !     S step
   call mat_zero(Uwrk); Swrk = 0.0_wp
   call apply_outerproduct(Uwrk, B, U)
   call mat_mult(Swrk, U1, Uwrk)         ! - Sdot
   !> Construct solution S1
   call mat_axpby(S, 1.0_wp, Swrk, -dt)


   !     L step
   call mat_zero(Uwrk)
   call mat_mult(Uwrk, U, transpose(S))  ! L0.T
   call apply_outerproduct(U, B, U1)     ! Ldot.T
   !> Construct solution Uwrk.T
   call mat_axpby(Uwrk, 1.0_wp, U, dt)
   !> Update S
   call mat_mult(S, Uwrk, U1)


   !> Copy data to output
   call mat_copy(U, U1)

   !> RKlib
   ! set initial condition
   Xtmp = reshape(X0, (/ nx**2, 1 /))
   call mat_zero(X_RKlib_mat)
   X_RKlib_mat(1)%state = Xtmp(:,1)
   ! initialize exponential propagator
   A_prop_mat = rklib_lyapunov_mat(dt)
   ! run step
   call A_prop_mat%matvec(X_RKlib_mat(1), X_RKlib_mat(2))
   ! recover output
   Xmat_rklib_out = reshape(X_RKlib_mat(2)%state, (/ nx, nx /))

   !> check result
   do i = 1, rk
      Xmat_DLRA_out(:,i) = U(i)%state
   end do
   USUT = matmul( Xmat_DLRA_out, matmul(S, transpose(Xmat_DLRA_out) ) )

   write(*,'(A26,F8.6,A32,F8.6)') ' L-T, single step dt = ', dt, &
            & '|| X_rk - X_DLRA ||_2 ~ ', norm_fro(Xmat_rklib_out - USUT)

   !------------------
   ! COMPUTE SINGLE M STEP + G STEP USING KEXPM FOR DIFFERENT DT AND COMPARE TO RKLIB
   !------------------

   do j = 1, 4
      dt = dt/10
   
      ! Reset input
      do i = 1,rk
         U(i)%state = Xmat_DLRA_in(:,i)
      end do
      S = S0

      ! M step
      torder = 1
      Ak_prop_vec = krylov_exptA_laplacian_vec(A,dt)
      call numerical_low_rank_splitting_integrator(U, S, Ak_prop_vec, B, dt, dt, torder, info)

      !> RKlib
      ! set initial condition
      Xtmp = reshape(X0, (/ nx**2, 1 /))
      call mat_zero(X_RKlib_mat)
      X_RKlib_mat(1)%state = Xtmp(:,1)
      ! initialize exponential propagator
      A_prop_mat = rklib_lyapunov_mat(dt)
      ! run step
      call A_prop_mat%matvec(X_RKlib_mat(1), X_RKlib_mat(2))
      ! recover output
      Xmat_rklib_out = reshape(X_RKlib_mat(2)%state, (/ nx, nx /))

      !> check result
      do i = 1, rk
         Xmat_DLRA_out(:,i) = U(i)%state
      end do
      USUT = matmul( Xmat_DLRA_out, matmul(S, transpose(Xmat_DLRA_out) ) )

      write(*,'(A26,F8.6,A32,F8.6)') ' L-T, single step dt = ', dt, &
            & '|| X_rk - X_DLRA ||_2 ~ ', norm_fro(Xmat_rklib_out - USUT)
   end do

   !------------------
   ! COMPUTE SECOND ORDER M STEP + G STEP USING KEXPM FOR DIFFERENT DT AND COMPARE TO RKLIB
   !------------------

   dt = 1.0_wp
   do j = 1, 5
      dt = dt/10
   
      ! Reset input
      do i = 1,rk
         U(i)%state = Xmat_DLRA_in(:,i)
      end do
      S = S0

      ! M step
      Ak_prop_vec = krylov_exptA_laplacian_vec(A,dt/2)
      torder = 2
      call numerical_low_rank_splitting_integrator(U, S, Ak_prop_vec, B, dt, dt, torder, info)

      !> RKlib
      ! set initial condition
      Xtmp = reshape(X0, (/ nx**2, 1 /))
      call mat_zero(X_RKlib_mat)
      X_RKlib_mat(1)%state = Xtmp(:,1)
      ! initialize exponential propagator
      A_prop_mat = rklib_lyapunov_mat(dt)
      ! run step
      call A_prop_mat%matvec(X_RKlib_mat(1), X_RKlib_mat(2))
      ! recover output
      Xmat_rklib_out = reshape(X_RKlib_mat(2)%state, (/ nx, nx /))

      !> check result
      do i = 1, rk
         Xmat_DLRA_out(:,i) = U(i)%state
      end do
      USUT = matmul( Xmat_DLRA_out, matmul(S, transpose(Xmat_DLRA_out) ) )

      write(*,'(A26,F8.6,A32,F8.6)') ' Strang, single step dt = ', dt, &
            & '|| X_rk - X_DLRA ||_2 ~ ', norm_fro(Xmat_rklib_out - USUT)
   end do

   !------------------
   ! COMPUTE DLRA FOR DIFFERENT TAU AND COMPARE TO RKLIB & STUART-BARTELS
   !------------------

   ! Reset input
   do i = 1,rk
      U(i)%state = Xmat_DLRA_in(:,i)
   end do
   S = S0

   Tend = 0.2_wp
   dt = 0.00001_wp

   ! DLRA
   Ak_prop_vec = krylov_exptA_laplacian_vec(A,dt)
   torder = 1

   do j = 1, 5
      ! run step
      call numerical_low_rank_splitting_integrator(U, S, Ak_prop_vec, B, Tend, dt, torder, info)

      !> check result
      do i = 1, rk
         Xmat_DLRA_out(:,i) = U(i)%state
      end do
      USUT = matmul( Xmat_DLRA_out, matmul(S, transpose(Xmat_DLRA_out) ) )
      call print_mat(nx, nx, USUT)

      write(*,'(A16,F8.4,A32,E16.8)') ' DRLA, Tend = ', j*Tend,&
            & '|| X_DLRA - X_ref ||_2 ~ ', norm_fro(Xmat - USUT)
   end do

   !> RKlib
   ! set initial condition
   Xtmp = reshape(X0, (/ nx**2, 1 /))
   call mat_zero(X_RKlib_mat)
   X_RKlib_mat(1)%state = Xtmp(:,1)

   ! initialize exponential propagator
   A_prop_mat = rklib_lyapunov_mat(Tend)

   Tend = 0.2_wp
   do j = 1, 5
      ! run step
      call A_prop_mat%matvec(X_RKlib_mat(1), X_RKlib_mat(2))
      ! recover output
      Xmat_rklib_out = reshape(X_RKlib_mat(2)%state, (/ nx, nx /))

      write(*,'(A16,F8.4,A32,E16.8)') ' RKlib, Tend = ', j*Tend,&
            & '|| X_rk   - X_ref ||_2 ~ ', norm_fro(Xmat - Xmat_rklib_out)
      
      ! continue computation
      X_RKlib_mat(1)%state = X_RKlib_mat(2)%state
   end do

   write(*,*) 'DONE'

end program demo