program demo
  use LightKrylov
  use Laplacian
  use Laplacian_Lyapunov
  use LightROM_LyapunovSolvers
  !use stdlib_io_npy, only : save_npy
  implicit none

  !------------------------------------------------
  !-----     LINEAR OPERATOR INVESTIGATED     -----
  !------------------------------------------------

  ! Laplacian in Lyapunov equation
  !> Exponential propagator.
  type(exponential_prop_lyap), allocatable :: A_lyap
  !> Sampling time.
  real(kind=wp), parameter :: tau = 0.1_wp
  !> Krylov subspace.
  type(state_vector_lyap), allocatable :: X_lyap(:)
  type(state_vector_lyap), allocatable :: wrk_lyap(:)

  ! Regular Laplacian
  !> Exponential propagator.
  type(exponential_prop),      allocatable :: A
  !> Number of eigenvalues we wish to converge.
  integer, parameter :: nev = 10
  !> Krylov subspace.
  type(state_vector),      allocatable :: X(:)
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
  integer, parameter :: torder = 1
  real(wp) :: Tend
  type(state_vector), allocatable :: B(:)
  type(state_vector), allocatable :: U(:)
  real(wp) :: S(rk,rk)
  real(wp) :: dt          !> timestep

  !> testing
  real(wp) :: Xmat_lyap(nx**2,2)

  call initialize_mesh

  A_lyap = exponential_prop_lyap(tau)
  A = exponential_prop(tau)

  allocate(X_lyap(1))
  allocate(wrk_lyap(1))

  !> Random initial Krylov vector.
  call random_number(X_lyap(1)%state); 
  !> symmetrize
  tmp = reshape(X_lyap(1)%state,shape(tmp)) + reshape(X_lyap(1)%state,shape(tmp),order=(/2,1/))
  X_lyap(1)%state = reshape(tmp,(/nx**2/))
  !> normalize
  alpha = X_lyap(1)%norm();
  call X_lyap(1)%scal(1.0_wp / alpha)

  Xmat_lyap(:,1) = X_lyap(1)%state
  !> Apply propagator
  call A_lyap%matvec(X_lyap(1),wrk_lyap(1))
  Xmat_lyap(:,2) = wrk_lyap(1)%state

  !write(*,*) 'Lypapunov equation:'
  !write(*,*) '    Initial condition:'
  !do i = 1,nx
  !    write(*,'(10F8.3)') Xmat_lyap((i-1)*nx+1:i*nx,1)
  !enddo
  !write(*,*) '    Output at T =',tau
  !do i = 1,nx
  !    write(*,'(10F8.3)') Xmat_lyap((i-1)*nx+1:i*nx,2)
  !enddo

  !> Initialize Krylov subspace.
  allocate(X(1:kdim+1)) ; call initialize_krylov_subspace(X)

  !> Random initial Krylov vector.
  call random_number(X(1)%state) ; alpha = X(1)%norm() ; call X(1)%scal(1.0_wp / alpha)

  !------------------------------------------
  !-----     EIGENVALUE COMPUTATION     -----
  !------------------------------------------

  !> Call to LightKrylov.
  call eigs(A, X, v, lambda, residuals, info, nev=nev)

  !> Transform eigenspectrum from unit-disk to standard complex plane.
  lambda = log(lambda) / tau

  !--------------------------------
  !-----     SAVE TO DISK     -----
  !--------------------------------

  !> Save the eigenspectrum.
  !call save_eigenspectrum(lambda%re, lambda%im, residuals, "example/laplacian1D/eigenspectrum.npy")
  write(*,*) 'Eigenspectrum Laplacian:'
  do i = 1,nev
      write(*,'(F8.3,A3,F8.3,A2)') lambda(i)%re, ' + ', lambda(i)%im,' i'
  enddo

  allocate(B(1))
  B(1)%state = 1.0
  Tend = 1.0
  !> Random initial condition of rank 4
  allocate(U(1:rk)); 
  do i = 1,rk
      call random_number(U(i)%state)
  enddo
  !> Orthonormalize (in-place)
  S = 0.0_wp
  call qr_factorization(U,S,info)
  dt = 0.1
  call numerical_low_rank_splitting_integrator(U, S, A, B, Tend, dt, torder, info)

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