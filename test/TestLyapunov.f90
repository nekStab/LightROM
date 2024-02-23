module TestLyapunov
   use LightKrylov
   use TestVector
   use TestMatrices
   Use LightROM_LyapunovUtils
   Use LightROM_wplib
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
          new_unittest("Dense Matrix Exponential", test_dense_matrix_exponential), &
          new_unittest("Development tests", playground) &
          ]
 
     return
   end subroutine collect_lyapunov_utils_testsuite

   subroutine test_dense_matrix_exponential(error)

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
      Eref = 0.0_wp
      forall (i=1:n) Eref(i, i) = 1.0_wp
      do i = 1, n-1
         do j = 1, n-i
            Eref(i,i+j) = m**j/factorial(j)
         end do
      end do
      
      E = 0.0_wp
      call expm(E, A)

      call check(error, maxval(E-Eref) < rtol)
      
      return
   end subroutine test_dense_matrix_exponential
   
   subroutine playground(error)

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
      integer, parameter :: nk = 20
      real(kind=wp) :: Xmat(test_size, nk), Qmat(test_size)
      real(kind=wp) :: Xrefmat(test_size)
      real(kind=wp) :: alpha
      real(kind=wp) :: Xreshape(test_size*kdim,1)
      real(kind=wp) :: Xmatr(kdim,test_size)
      real(wp) :: pad(1)
      real(wp) :: tau
      real(wp) :: difference(nk)

      pad = 0.0_wp

      ! --> Initialize matrix.
      A = rmatrix() ; call random_number(A%data)
      !A%data = reshape((/  1.0, 0.0, 0.0, &
      !                    &0.0, 1.0, 0.0, &
      !                    &0.0, 0.0, 1.0 /), (/ kdim, kdim /) )

      allocate(Q(1)) ; call random_number(Q(1)%data)
      Qmat = Q(1)%data
      allocate(Xref(1)) ; call mat_zero(Xref)
      allocate(Xkryl(1:nk)) ; call mat_zero(Xkryl)

      Amat = 0.0_wp;
      !forall (i=1:kdim) Amat(i, i) = 1.0_wp
      Amat = A%data
      !do i = 1, kdim
      !   do j = 1, kdim
      !      call random_number(A(i,j))
      !   end do 
      !end do

      tau = 0.1_wp

      !write(*,*) 'Amat test'
      !do i = 1, kdim
      !   write(*,'(10E15.8)') Amat(i,1:kdim)
      !enddo
      write(*,*) 'Qmat test'
      do i = 1, kdim
         write(*,'(1E15.8)') Qmat(i)
      enddo

      Emat = 0.0_wp
      Xrefmat = 0.0_wp
      call expm(Emat, tau*Amat)
      do i = 1, test_size
         do j = 1,kdim
            Xrefmat(i) = Xrefmat(i) + Emat(i,j) * Qmat(j)
         end do
      end do
      Xref(1)%data = Xrefmat

      !write(*,*) 'Emat test'
      !do i = 1, kdim
      !   write(*,'(10E15.8)') Emat(i,1:kdim)
      !enddo

      Emat = 0.0_wp
      do i = 1, nk
         call kexpm(Xkryl(i:i), A, Q(1:1), tau, 1.0e-6_wp, info, nkryl = i)
         call Xkryl(i)%axpby(1.0_wp, Xref(1), -1.0_wp)
         difference(i) = Xkryl(i)%norm()
      end do


      do i = 1,nk
         write(*,'(a13,2E15.8)') '|| dx ||_2 = ', difference(i)
      enddo

      STOP 1

      !write(*,*) A%data

      ! --> Initialize Krylov subspace.
      allocate(Q(1:kdim)) ; call random_number(Q(1)%data)
      alpha = Q(1)%norm() ; call Q(1)%scal(1.0D+00 / alpha)
      do k = 2, size(Q)
!         call Q(k)%zero()
         call random_number(Q(k)%data)
      enddo

      ! Copy data
      !do k = 1, kdim
      !   Qmat(:, k) = Q(k)%data
      !enddo

      !Amat = matmul(A, Qmat)
      allocate(Xref(1:kdim));
      call mat_mult(Xref,Q,Amat)
      ! Copy data
      do k = 1, kdim
         Xmat(:, k) = Xref(k)%data
      enddo

      !write(*,*) 'A'
      !write(*,*) shape(A)
      !write(*,*) 'A'
      !do i = 1, kdim
      !   write(*,'(3F8.3)') A(i,1:kdim)
      !enddo
      !do i = 1, test_size
      !   write(*,'(10F8.3)') Xmat(i,1:kdim)
      !enddo
      !write(*,*) ''
      !Xreshape = reshape(Xmat,shape(Xreshape))
      !!write(*,'(20F8.3)') Xreshape(1:test_size*kdim,1)
      !Xmatr = reshape(Xmat,(/kdim, test_size/), pad, (/ 2,1 /))
      !do i = 1, kdim
      !   write(*,'(10F8.3)') Xmatr(i,1:test_size)
      !enddo
      !Xreshape = reshape(Xmat,(/test_size*kdim, 1/), pad, (/ 2,1 /))
      !write(*,'(20F8.3)') Xreshape(1:test_size*kdim,1)
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