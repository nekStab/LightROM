module TestLyapunov
   use LightKrylov
   use TestVector
   use TestMatrices
   Use LightROM_LyapunovUtils
   use testdrive  , only : new_unittest, unittest_type, error_type, check
   use stdlib_math, only : all_close
   implicit none
 
   private
 
   public :: collect_lyapunov_utils_testsuite

  contains
 
   !-------------------------------------------
   !-----                                 -----
   !-----     TEST SUITE FOR THE DLRA     -----
   !-----                                 -----
   !-------------------------------------------
 
   subroutine collect_lyapunov_utils_testsuite(testsuite)
     !> Collection of tests.
     type(unittest_type), allocatable, intent(out) :: testsuite(:)
 
     testsuite = [&
          new_unittest("Matrix product direct", test_direct_krylov_matrix_product), &
          new_unittest("Matrix product transpose", test_transpose_krylov_matrix_product), &
          new_unittest("Development tests", playground) &
          ]
 
     return
   end subroutine collect_lyapunov_utils_testsuite

   subroutine test_direct_krylov_matrix_product(error)
      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      !> Test bases.
      class(rvector), dimension(:), allocatable :: A(:)
      class(rvector), dimension(:), allocatable :: C(:)
      !> Krylov subspace dimension.
      integer, parameter :: kdim1 = 3
      !> Number of columns in coefficien matrix
      integer, parameter :: kdim2 = 4
      !> Test matrices.
      real(kind=wp)               :: B(kdim1, kdim2)
      real(kind=wp)               :: Amat(test_size, kdim1)
      real(kind=wp)               :: Cmat(test_size, kdim2)
      !> Misc.
      integer :: i,j

      !> Initialize basis and copy data to matrix
      allocate(A(1:kdim1))
      do i = 1, size(A)
         call random_number(A(i)%data)
         Amat(:, i) = A(i)%data
      enddo
      allocate(C(1:kdim2))
      do i = 1, size(C)
         call C(i)%zero()
      enddo
      B = 0.0_wp
      do i = 1, size(A)
         do j = 1, size(C)
            call random_number(B(i,j))
         enddo
      enddo
      !> Compute product
      call mat_mult(C,A,B)
      !> Copy data
      do i = 1, kdim2
         Cmat(:, i) = C(i)%data
      enddo
      call check(error, all_close(matmul(Amat, B), Cmat, rtol, atol) )
      return
   end subroutine test_direct_krylov_matrix_product

   subroutine test_transpose_krylov_matrix_product(error)
      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      !> Test bases.
      class(rvector), dimension(:), allocatable :: A(:)
      class(rvector), dimension(:), allocatable :: B(:)
      !> Krylov subspace dimension.
      integer, parameter :: kdim1 = 3
      integer, parameter :: kdim2 = 4
      !> Test matrices.
      real(kind=wp)               :: C(kdim1, kdim2)
      real(kind=wp)               :: Amat(test_size, kdim1)
      real(kind=wp)               :: Bmat(test_size, kdim2)
      !> Misc.
      integer :: k

      !> Initialize bases and copy data to matrices
      allocate(A(1:kdim1))
      Amat = 0.0_wp
      do k = 1, size(A)
         call random_number(A(k)%data)
         Amat(:, k) = A(k)%data
      enddo
      allocate(B(1:kdim2))
      Bmat = 0.0_wp
      do k = 1, size(B)
         call random_number(B(k)%data)
         Bmat(:, k) = B(k)%data
      enddo
      C = 0.0_wp
      !> Compute product
      call mat_mult(C,A,B)
      call check(error, all_close(matmul(transpose(Amat), Bmat), C, rtol, atol) )
      return
   end subroutine test_transpose_krylov_matrix_product
   
   subroutine playground(error)

      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      !> Basis vectors.
      class(rvector), allocatable :: Q(:)
      class(rvector), allocatable :: X(:)
      !> Krylov subspace dimension.
      integer, parameter :: kdim = 3
      !> Test matrix.
      real(kind=wp) :: A(kdim, kdim)
      !> GS factors.
      real(kind=wp) :: R(kdim, kdim)
      real(kind=wp) :: Id(kdim, kdim)
      !> Information flag.
      integer :: info
      !> Misc.
      integer :: i,j,k
      real(kind=wp) :: Xmat(test_size, kdim), Qmat(test_size, kdim)
      real(kind=wp) :: Dmat(test_size, kdim)
      real(kind=wp) :: alpha

      ! --> Initialize matrix.
      !A = rmatrix() ; call random_number(A%data)
      do i = 1,kdim
         do j = 1,kdim
            call random_number(A(i,j))
         enddo
      enddo
      !write(*,*) 'A'
      !do i = 1, kdim
      !   write(*,'(3F8.3)') A(i,1:kdim)
      !enddo

      !write(*,*) A%data

      ! --> Initialize Krylov subspace.
      allocate(Q(1:kdim)) ; call random_number(Q(1)%data)
      alpha = Q(1)%norm() ; call Q(1)%scal(1.0D+00 / alpha)
      do k = 2, size(Q)
!         call Q(k)%zero()
         call random_number(Q(k)%data)
      enddo

      ! Copy data
      do k = 1, kdim
         Qmat(:, k) = Q(k)%data
      enddo

      !Amat = matmul(A, Qmat)
      allocate(X(1:kdim));
      call mat_mult(X,Q,A)
      ! Copy data
      do k = 1, kdim
         Xmat(:, k) = X(k)%data
      enddo

      !write(*,*) 'A'
      !do i = 1, kdim
      !   write(*,'(3F8.3)') A(i,1:kdim)
      !enddo
      !write(*,*) 'Q'
      !do i = 1, test_size
      !   write(*,'(3F8.3)') Qmat(i,1:kdim)
      !enddo
      !write(*,*) 'X = Q @ A'
      !do i = 1, test_size
      !   write(*,'(3F8.3)') Xmat(i,1:kdim)
      !enddo
      !Dmat = matmul(Qmat,A)
      !write(*,*) 'X = Q @ A'
      !do i = 1, test_size
      !   write(*,'(3F8.3)') Dmat(i,1:kdim)
      !enddo

      !write(*,*) size(A) ! not possible!
      !write(*,*) 'Size Q (n x r) abstract vector'
      !write(*,*) size(Q)
      !write(*,*) 'Size Xmat (n x r)'
      !write(*,*) size(Xmat), size(Xmat,1), size(Xmat, 2)
      
      call check(error, 0.0_wp < rtol)
      
      return
   end subroutine playground

end module TestLyapunov