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
          new_unittest("Krylov Matrix Exponential", test_krylov_matrix_exponential) &
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
      !> Misc.
      integer :: i,j,k
      integer, parameter :: nk = 10
      real(kind=wp) :: Xmat(test_size), Qmat(test_size)
      real(kind=wp) :: alpha
      real(wp) :: tau

      ! --> Initialize matrix.
      A = rmatrix() ; call random_number(A%data)
      
      allocate(Q(1)) ; call random_number(Q(1)%data)
      Qmat = Q(1)%data
      allocate(Xref(1)) ; call mat_zero(Xref)
      allocate(Xkryl(1)) ; call mat_zero(Xkryl)

      Amat = 0.0_wp;
      Amat = A%data

      tau = 0.1_wp

      Emat = 0.0_wp
      Xmat = 0.0_wp
      !> Comparison is dense computation (10th order Pade approximation)
      call expm(Emat, tau*Amat)
      Xmat = matmul(Emat,Qmat)
      !> Copy reference data into Krylov vector
      Xref(1)%data = Xmat

      !> Compute Krylov matrix exponential for different krylov subspace sizes
      do i = 1, nk
         call kexpm(Xkryl(1:1), A, Q(1:1), tau, 1.0e-6_wp, info, nkryl = i)
         !> Compute 2-norm of the error
         call Xkryl(1)%axpby(1.0_wp, Xref(1), -1.0_wp)
         alpha = Xkryl(1)%norm()
         write(*,'(a12,i2,a12,E15.7)') 'Krylov step ', i, ': error_2 = ', alpha
      end do
     
      call check(error, alpha < rtol)
      
      return
   end subroutine test_krylov_matrix_exponential

end module TestExpm