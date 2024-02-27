module TestExpm
   use LightKrylov
   use TestVector
   use TestMatrices
   Use LightROM_expmlib
   use testdrive  , only : new_unittest, unittest_type, error_type, check
   use stdlib_math, only : all_close
   implicit none
 
   private
 
   public :: collect_expm_testsuite

  contains
 
   !---------------------------------------------------------
   !-----                                               -----
   !-----     TEST SUITE FOR THE MATRIX EXPONENTIAL     -----
   !-----                                               -----
   !---------------------------------------------------------
 
   subroutine collect_expm_testsuite(testsuite)
     !> Collection of tests.
     type(unittest_type), allocatable, intent(out) :: testsuite(:)
 
     testsuite = [&
          new_unittest("Dense Matrix Exponential", test_dense_matrix_exponential), &
          new_unittest("Krylov Matrix Exponential", test_krylov_matrix_exponential), &
          new_unittest("Block Krylov Matrix Exponential", test_block_krylov_matrix_exponential) &
          ]
 
     return
   end subroutine collect_expm_testsuite

   subroutine test_dense_matrix_exponential(error)
      !> This function tests the scaling and squaring followed by rational Pade approximation
      ! of the matrix exponential for a matrix for which the exponential propagator is known
      ! analytically

      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      !> Problem dimension.
      integer, parameter :: n = 5
      integer, parameter :: m = 6
      !> Test matrix.
      real(kind=wp) :: A(n, n)
      real(kind=wp) :: E(n, n)
      real(kind=wp) :: Eref(n, n)
      integer :: i, j

      ! --> Initialize matrix.
      A = 0.0_wp
      do i = 1, n-1
         A(i,i+1) = m*1.0_wp
      end do
      ! --> Reference with analytical exponential
      Eref = 0.0_wp
      forall (i=1:n) Eref(i, i) = 1.0_wp
      do i = 1, n-1
         do j = 1, n-i
            Eref(i,i+j) = Eref(i,i+j-1)*m/j
         end do
      end do
      ! --> Compute exponential numerically
      E = 0.0_wp
      call expm(E, A)

      call check(error, maxval(E-Eref) < rtol)
      
      return
   end subroutine test_dense_matrix_exponential

   subroutine test_krylov_matrix_exponential(error)
      !> This function tests the Krylov based approximation of the action of the exponential
      ! propagator against the dense computation for a random operator, a random RHS and a 
      ! typical value of tau.

      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      class(rmatrix), allocatable :: A
      !> Basis vectors.
      class(rvector), allocatable :: Q(:)
      class(rvector), allocatable :: Xref(:)
      class(rvector), allocatable :: Xkryl(:)
      !> Krylov subspace dimension.
      integer, parameter :: kdim = test_size
      !> Test matrix.
      real(kind=wp) :: Amat(kdim, kdim)
      real(kind=wp) :: Emat(kdim, kdim)
      !> GS factors.
      real(kind=wp) :: R(kdim, kdim)
      real(kind=wp) :: Id(kdim, kdim)
      !> Information flag.
      integer :: info
      !> Test parameters
      integer, parameter         :: nkmax = 15
      integer, parameter         :: p     = 1
      real(kind=wp), parameter   :: tau   = 0.1_wp
      real(kind=wp), parameter   :: tol   = 1e-10_wp
      !> Misc.
      integer :: i,j,k
      real(kind=wp) :: Xmat(test_size,p), Qmat(test_size,p)
      real(kind=wp) :: alpha
      real(kind=wp) :: err(p,p)

      Amat = 0.0_wp; Emat = 0.0_wp; Xmat = 0.0_wp
      allocate(Xref(1:p)) ; call mat_zero(Xref)
      allocate(Xkryl(1:p)) ; call mat_zero(Xkryl)

      ! --> Initialize operator.
      A = rmatrix() ; call random_number(A%data)
      Amat = A%data     
      ! --> Initialize rhs.
      allocate(Q(1:p)) ; 
      do i = 1,p
         call random_number(Q(i)%data)
         Qmat(:,i) = Q(i)%data
      end do
      
      !> Comparison is dense computation (10th order Pade approximation)
      call expm(Emat, tau*Amat)
      Xmat = matmul(Emat,Qmat)
      !> Copy reference data into Krylov vector
      do i = 1,p
         Xref(i)%data = Xmat(:,i)
      end do

      !> Compute Krylov matrix exponential for different krylov subspace sizes
      call kexpm(Xkryl(1:p), A, Q(1:p), tau, tol, info, verbosity = .true., nkryl = nkmax)
      do i = 1,p
         call Xkryl(i)%axpby(1.0_wp, Xref(i), -1.0_wp)
      end do
      !> Compute 2-norm of the error
      call mat_mult(err,Xkryl(1:p),Xkryl(1:p))
      alpha = sqrt(norm_fro(err))
      write(*, *) '    true error:          ||error||_2 = ', alpha
     
      call check(error, alpha < rtol)
      
      return
   end subroutine test_krylov_matrix_exponential

   subroutine test_block_krylov_matrix_exponential(error)
      !> This function tests the Krylov based approximation of the action of the exponential
      ! propagator against the dense computation for a random operator, a random RHS and a 
      ! typical value of tau.

      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      class(rmatrix), allocatable :: A
      !> Basis vectors.
      class(rvector), allocatable :: Q(:)
      class(rvector), allocatable :: Xref(:)
      class(rvector), allocatable :: Xkryl(:)
      !> Krylov subspace dimension.
      integer, parameter :: kdim = test_size
      !> Test matrix.
      real(kind=wp) :: Amat(kdim, kdim)
      real(kind=wp) :: Emat(kdim, kdim)
      !> GS factors.
      real(kind=wp) :: R(kdim, kdim)
      real(kind=wp) :: Id(kdim, kdim)
      !> Information flag.
      integer :: info
      !> Test parameters
      integer, parameter         :: nkmax = 15
      integer, parameter         :: p     = 2
      real(kind=wp), parameter   :: tau   = 0.1_wp
      real(kind=wp), parameter   :: tol   = 1e-10_wp
      !> Misc.
      integer :: i,j,k
      real(kind=wp) :: Xmat(test_size,p), Qmat(test_size,p)
      real(kind=wp) :: alpha
      real(kind=wp) :: err(p,p)

      Amat = 0.0_wp; Emat = 0.0_wp; Xmat = 0.0_wp
      allocate(Xref(1:p)) ; call mat_zero(Xref)
      allocate(Xkryl(1:p)) ; call mat_zero(Xkryl)

      ! --> Initialize operator.
      A = rmatrix() ; call random_number(A%data)
      Amat = A%data     
      ! --> Initialize rhs.
      allocate(Q(1:p)) ; 
      do i = 1,p
         call random_number(Q(i)%data)
         Qmat(:,i) = Q(i)%data
      end do
      
      !> Comparison is dense computation (10th order Pade approximation)
      call expm(Emat, tau*Amat)
      Xmat = matmul(Emat,Qmat)
      !> Copy reference data into Krylov vector
      do i = 1,p
         Xref(i)%data = Xmat(:,i)
      end do

      !> Compute Krylov matrix exponential for different krylov subspace sizes
      call kexpm(Xkryl(1:p), A, Q(1:p), tau, tol, info, verbosity = .true., nkryl = nkmax)
      do i = 1,p
         call Xkryl(i)%axpby(1.0_wp, Xref(i), -1.0_wp)
      end do
      !> Compute 2-norm of the error
      call mat_mult(err,Xkryl(1:p),Xkryl(1:p))
      alpha = sqrt(norm_fro(err))
      write(*, *) '    true error:          ||error||_2 = ', alpha
     
      call check(error, alpha < rtol)
      
      return
   end subroutine test_block_krylov_matrix_exponential

end module TestExpm