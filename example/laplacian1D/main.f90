program demo
  use LightKrylov
  use Laplacian
  use Laplacian_Lyapunov
  use LightROM_LyapunovSolvers
  use LightROM_utils
  use LightROM_expmlib
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
  real(kind=wp), parameter :: tau = 0.01_wp !0.1_wp
  real(kind=wp), parameter :: tau_long = 0.1_wp !1.0_wp

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
  integer, parameter :: rk = 4
  integer, parameter :: torder = 2
  real(wp) :: Tend
  real(wp) :: S(rk,rk)
  real(wp) :: dt          !> timestep

  !> testing laplace operator
  real(wp)  :: Xvec(nx,4)
  real(wp)  :: Avec(nx,nx)
  real(wp)  :: Evec(nx,nx)
  !> testing Lyapunov equation
  real(wp)  :: Xmat(nx**2,5)
  real(wp)  :: rXmat(nx,nx)
  real(wp)  :: A1mat(nx**2,nx**2)
  real(wp)  :: A2mat(nx**2,nx**2)
  real(wp)  :: Amat(nx**2,nx**2)
  real(wp)  :: Emat(nx**2,nx**2)
  real(wp)  :: q_flat(nx**2)
  real(wp)  :: Xtmp(nx**2,1)
  !> testing DLRA
  real(wp)  :: Xmat_lyap(nx**2,3)
  real(wp)  :: Xmat_DLRA_0(nx,rk)
  real(wp)  :: Xmat_DLRA_in(nx,rk)
  real(wp)  :: Xmat_DLRA_out(nx,rk)
  type(state_vector), allocatable :: B(:)
  type(state_vector), allocatable :: U(:)
  type(state_vector), allocatable :: Uwrk(:)
  real(wp)  :: US(nx,rk)
  real(wp)  :: USUT(nx,nx)
  type(state_matrix), allocatable :: X_DLRA_mat(:)
  real(wp)  :: Xmat_rklib_out(nx, rk)

  ! Laplace operator
  call initialize_mesh

  Avec = 0.0_wp
  forall (i=1:nx)   Avec(i,i)   = -2.0/dx**2
  forall (i=1:nx-1) Avec(i+1,i) = 1.0/dx**2
  forall (i=1:nx-1) Avec(i,i+1) = 1.0/dx**2
  
  !write(*,*) 'Avec - Laplace Operator'
  !call print_mat(nx, nx, Avec)

  Xvec = 0.0_wp
  allocate(X_vec(1:4))
  !> Random initial Krylov vector.
  call random_number(X_vec(1)%state)
  !X_vec(1)%state = (/ 1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp /)
  Xvec(:,1) = X_vec(1)%state

  Xvec(:,2) = matmul(Avec,Xvec(:,1))
  call A%matvec(X_vec(1), X_vec(2))
  Xvec(:,3) = X_vec(2)%state

  !write(*,*) 'Xvec - Laplace'
  !call print_mat(nx, 3, Xvec)
  
  !> Initialize exponential propagators for the laplace equation.
  A_prop_vec  = rklib_exptA_laplacian_vec(tau)
  Ak_prop_vec = krylov_exptA_laplacian_vec(A,tau)

  Evec = 0.0_wp
  call expm(Evec,tau*Avec)
  !write(*,*) 'Evec - expm(tau*Avec) Laplace Operator'
  !call print_mat(nx, nx, Evec)

  call X_vec(2)%zero()
  call X_vec(3)%zero()
  call A_prop_vec%matvec(X_vec(1), X_vec(2))
  call Ak_prop_vec%matvec(X_vec(1), X_vec(3))
  X_vec(4)%state = matmul(Evec,Xvec(:,1))
  do i = 2,4
   Xvec(:,i) = X_vec(i)%state
  end do

  !write(*,*) '-----------------------------------------------------'
  !write(*,*) '       expm(tau*Avec) @ X(0), Laplace Operator'
  !write(*,*) '               tau = ', tau
  !write(*,*) '-----------------------------------------------------'
  !write(*,*) '    X(0)         RKlib         Krylov        Pade'
  !write(*,*) '-----------------------------------------------------'
  !call print_mat(nx, 4, Xvec)
  call X_vec(2)%axpby(1.0_wp, X_vec(4), -1.0_wp)
  call X_vec(3)%axpby(1.0_wp, X_vec(4), -1.0_wp)
  !write(*,"(A18,3(E12.6,'  '))") '|error|_2     ',X_vec(2)%norm(),X_vec(3)%norm(),0.0

  call expm(Evec,tau_long*Avec)
  A_prop_vec  = rklib_exptA_laplacian_vec(tau_long)
  Ak_prop_vec = krylov_exptA_laplacian_vec(A,tau_long)
  call X_vec(2)%zero()
  call X_vec(3)%zero()
  call A_prop_vec%matvec(X_vec(1), X_vec(2))
  call Ak_prop_vec%matvec(X_vec(1), X_vec(3))
  X_vec(4)%state = matmul(Evec,Xvec(:,1))
  do i = 2,4
   Xvec(:,i) = X_vec(i)%state
  end do
  !write(*,*) '-----------------------------------------------------'
  !write(*,*) '               tau = ', tau_long
  !write(*,*) '-----------------------------------------------------'
  !write(*,*) '    X(0)         RKlib         Krylov        Pade'
  !write(*,*) '-----------------------------------------------------'
  !call print_mat(nx, 4, Xvec)
  call X_vec(2)%axpby(1.0_wp, X_vec(4), -1.0_wp)
  call X_vec(3)%axpby(1.0_wp, X_vec(4), -1.0_wp)
  !write(*,"(A18,3(E12.6,'  '))") '|error|_2     ',X_vec(2)%norm(),X_vec(3)%norm(),0.0
  !write(*,*) '-----------------------------------------------------'
  !write(*,*)

  call initialize_mesh_matrix_eq

  ! A @ X   = kron(A,I) @ vec(X)
  A1mat = 0.0_wp
  forall (i=1:nx)    A1mat((i-1)*nx+1:i*nx,(i-1)*nx+1:i*nx) = Avec(1:nx,1:nx)
  ! (A @ X.T).T = kron(I,A) @ vec(X)
  A2mat = 0.0_wp
  forall (i=1:nx**2)     A2mat(i,i)    = -2.0/dx**2
  forall (i=1:nx*(nx-1)) A2mat(nx+i,i) =  1.0/dx**2
  forall (i=1:nx*(nx-1)) A2mat(i,nx+i) =  1.0/dx**2
  ! ( Avec @ X + X @ Avec.T ) = Amat @ vec(X)
  Amat = A1mat + A2mat
  
  !call print_mat(nx**2, nx**2, A1mat)
  !call print_mat(nx**2, nx**2, A2mat)
  !call print_mat(nx**2, nx**2, A2Tmat)
  !call print_mat(nx**2, nx**2, Amat)

  Xmat = 0.0_wp
  q_flat = 1.0_wp
  
  allocate(X_mat(1:4))
  call random_number(X_mat(1)%state)
  Xmat(:,1) = X_mat(1)%state
  rXmat = reshape(Xmat(:,1), (/ nx,nx /))
  !Xmat(:,2) = matmul(A1mat,Xmat(:,1)); Xmat(:,3) = reshape(matmul(Amat,rXmat), shape(Xmat(:,1)))
  !Xmat(:,2) = matmul(A2mat,Xmat(:,1)); Xmat(:,3) = reshape(matmul(rXmat,Amat), shape(Xmat(:,1)))
  Xmat(:,2) = matmul(Amat,Xmat(:,1)) + q_flat
  Xmat(:,3) = reshape(matmul(Avec,rXmat) + matmul(rXmat,Avec), shape(Xmat(:,1))) + q_flat

  call X_mat(2)%zero()
  call X_vec(3)%zero()
  call A_lyap%matvec(X_mat(1), X_mat(2))
  Xmat(:,4) = X_mat(2)%state

  !write(*,*) 'Xmat - Lyapunov equation for Laplace operator'
  !call print_mat(nx**2, 4, Xmat)

  !> Initialize exponential propagators for the Lyapunov equation of the Laplace operator
  A_prop_mat  = rklib_lyapunov_mat(tau)

  call X_mat(2)%zero()
  call X_mat(3)%zero()
  call A_prop_mat%matvec(X_mat(1), X_mat(2))
  do i = 2,4
   Xmat(:,i) = X_mat(i)%state
  end do
  
  A_prop_mat  = rklib_lyapunov_mat(tau_long)

  call X_mat(2)%zero()
  call X_mat(3)%zero()
  call A_prop_mat%matvec(X_mat(1), X_mat(2))
  do i = 2,4
   Xmat(:,i) = X_mat(i)%state
  end do
  
!
  !!> symmetrize
  !tmp = reshape(X_lyap(1)%state,shape(tmp)) + reshape(X_lyap(1)%state,shape(tmp),order=(/2,1/))
  !X_lyap(1)%state = reshape(tmp,(/nx**2/))
  !!> normalize
  !alpha = X_lyap(1)%norm();
  !call X_lyap(1)%scal(1.0_wp / alpha)
!
  !Xmat_lyap(:,1) = X_lyap(1)%state
  !!> Apply propagator using RKlib
  !call A_lyap%matvec(X_lyap(1),wrk_lyap(1))
  !Xmat_lyap(:,2) = wrk_lyap(1)%state
  !!> Apply propagator using kexpm
  !call A_kryl_lyap%matvec(X_lyap(1),wrk_lyap(1))
  !Xmat_lyap(:,3) = wrk_lyap(1)%state
!
  !write(*,*) 'Lypapunov equation:'
  !write(*,*) '    Initial condition:'
  !do i = 1,nx
  !    write(*,'(10F8.3)') Xmat_lyap((i-1)*nx+1:i*nx,1)
  !enddo
  !write(*,*) '    Output at T =',tau
  !do i = 1,nx
  !    write(*,'(10F8.3)') Xmat_lyap((i-1)*nx+1:i*nx,2)
  !enddo
  !write(*,*) '    Output at T =',tau
  !do i = 1,nx
  !    write(*,'(10F8.3)') Xmat_lyap((i-1)*nx+1:i*nx,3)
  !enddo
  !
  !!> Initialize Krylov subspace.
  !!allocate(X(1:kdim+1)) ; call initialize_krylov_subspace(X)

  !> Random initial Krylov vector.
  !call random_number(X(1)%state) ; alpha = X(1)%norm() ; call X(1)%scal(1.0_wp / alpha)

  !------------------------------------------
  !-----     EIGENVALUE COMPUTATION     -----
  !!------------------------------------------
!
  !!> Call to LightKrylov.
  !call eigs(A, X, v, lambda, residuals, info, nev=nev)
!
  !!> Transform eigenspectrum from unit-disk to standard complex plane.
  !lambda = log(lambda) / tau
!
  !!--------------------------------
  !!-----     SAVE TO DISK     -----
  !!--------------------------------
!
  !!> Save the eigenspectrum.
  !!call save_eigenspectrum(lambda%re, lambda%im, residuals, "example/laplacian1D/eigenspectrum.npy")
  !write(*,*) 'Eigenspectrum Laplacian:'
  !do i = 1,nev
  !    write(*,'(F8.3,A3,F8.3,A2)') lambda(i)%re, ' + ', lambda(i)%im,' i'
  !enddo
!
  allocate(B(1))
  B(1)%state = 0.0_wp

  Tend = 0.01_wp; dt = 0.01_wp
  Ak_prop_vec = krylov_exptA_laplacian_vec(A,dt/2)
  
  !> Random initial condition of rank 4
  Xmat_DLRA_in = 0.0_wp
  allocate(U(1:rk)); allocate(Uwrk(1:rk)); call mat_zero(Uwrk)
  U(1)%state = (/ 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp /)
  U(2)%state = (/ 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp /)
  U(3)%state = (/ 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp /)
  U(4)%state = (/ 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp /)
  do i = 1,rk
      !call random_number(U(i)%state)
      Xmat_DLRA_0(:,i) = U(i)%state
  enddo
  !> Orthonormalize (in-place)
  S = 0.0_wp
  call qr_factorization(U,S,info)
  do i = 1,rk
      Xmat_DLRA_in(:,i) = U(i)%state
  end do
  write(*,*) '    Input'
  call print_mat(nx, rk, Xmat_DLRA_0)
  write(*,*) '    U'
  call print_mat(nx, rk, Xmat_DLRA_in)
  write(*,*) '    S'
  call print_mat(rk, rk, S)
  call numerical_low_rank_splitting_integrator(U, S, Ak_prop_vec, B, Tend, dt, 2, info)
  write(*,*) '    S after DLRA'
  call print_mat(rk, rk, S)
  do i = 1,rk
     Xmat_DLRA_out(:,i) = U(i)%state
  end do
  write(*,*) '    U at T =',Tend
  call print_mat(nx, rk, Xmat_DLRA_out)
  US = matmul(Xmat_DLRA_out, S)
  USUT = 0.0_wp
  USUT = matmul(US, transpose(Xmat_DLRA_out))
  
  write(*,*) '    DLRA output at T =',Tend
  !call print_mat(nx, rk, Xmat_DLRA_out)
  call print_mat(nx, nx, USUT)

  ! RKlib
  allocate(X_DLRA_mat(1:2))
  Xtmp = reshape(Xmat_DLRA_0, (/ nx**2, 1 /))
  X_DLRA_mat(1)%state = Xtmp(:,1)
  !write(*,*) reshape(X_DLRA_mat(1)%state, (/ nx, rk /))

  ! Initialize exponential propagator
  A_prop_mat  = rklib_lyapunov_mat(dt)
  call A_prop_mat%matvec(X_DLRA_mat(1), X_DLRA_mat(2))
  Xmat_rklib_out = reshape(X_DLRA_mat(2)%state, (/ nx, rk /))
  write(*,*) '    RKlib output at T =',Tend
  call print_mat(nx, rk, Xmat_rklib_out)

  Xmat_rklib_out = Xmat_rklib_out - USUT

  write(*,*) norm_fro(Xmat_rklib_out)

  !> Reconstruct the leading eigenvectors from the Krylov basis.
  !do i = 1, nev
  !   !> Real part.
  !   call get_vec(wrk, X(1:kdim), v(:, i)%re)
  !   select type(wrk)
  !   type is(state_vector)
  !      eigenvectors(:, i)%re = wrk%state(1:nx)
  !   end select
!
  !   !> Imaginary part.
  !   call get_vec(wrk, X(1:kdim), v(:, i)%im)
  !   select type(wrk)
  !   type is(state_vector)
  !      eigenvectors(:, i)%im = wrk%state(1:nx)
  !   end select
  !enddo

  !> Save eigenvectors to disk.
  !call save_npy("example/laplacian1D/eigenvectors.npy", eigenvectors)
  !write(*,*) 'Eigenvectors (real part):'
  !do i = 1,nx
  !    write(*,'(10F8.3)') Xmat_lyap((i-1)*nx+1:i*nx,2)
  !enddo

  write(*,*) 'DONE'

end program demo