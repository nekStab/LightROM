module TestLyapunov
   use LightKrylov
   use TestVector
   use TestMatrices
   Use LightROM_LyapunovUtils
   Use LightROM_expmlib
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
          new_unittest("Krylov Matrix Exponential", test_krylov_matrix_exponential) &!, &
          ! new_unittest("Development tests", playground) &
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
   
   subroutine playground(error)

      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      class(rmatrix), allocatable :: A
      !> Basis vectors.
      class(rvector), allocatable :: Q(:)
      class(rvector), allocatable :: Xref(:)
      class(rvector), allocatable :: Xkryl(:)
      class(rvector), allocatable :: Xkrylc(:)
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
      real(kind=wp) :: Xmat(test_size, nk), Qmat(test_size)
      real(kind=wp) :: Xrefmat(test_size)
      real(kind=wp) :: alpha
      real(kind=wp) :: Xreshape(test_size*kdim,1)
      real(kind=wp) :: Xmatr(kdim,test_size)
      real(wp) :: pad(1)
      real(wp) :: tau, z, c
      real(wp) :: difference(nk,2)

      pad = 0.0_wp

      ! --> Initialize matrix.
      A = rmatrix() ; call random_number(A%data)
      !A%data = reshape((/  1.0, 0.0, 0.0, &
      !                    &0.0, 1.0, 0.0, &
      !                    &0.0, 0.0, 1.0 /), (/ kdim, kdim /) )

      !z = 0.0_wp
      !c = 0.5_wp
      !A%data = reshape((/ 0.04761905_wp,  c, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, &
      !                  & -c, 0.04761905_wp, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, &
      !                  & z, z, 0.14285714_wp,  c, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, &
      !                  & z, z, -c, 0.14285714_wp, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, &
      !                  & z, z, z, z, 0.23809524_wp,  c, z, z, z, z, z, z, z, z, z, z, z, z, z, z, &
      !                  & z, z, z, z, -c, 0.23809524_wp, z, z, z, z, z, z, z, z, z, z, z, z, z, z, &
      !                  & z, z, z, z, z, z, 0.33333333_wp,  c, z, z, z, z, z, z, z, z, z, z, z, z, &
      !                  & z, z, z, z, z, z, -c, 0.33333333_wp, z, z, z, z, z, z, z, z, z, z, z, z, &
      !                  & z, z, z, z, z, z, z, z, 0.42857143_wp,  c, z, z, z, z, z, z, z, z, z, z, &
      !                  & z, z, z, z, z, z, z, z, -c, 0.42857143_wp, z, z, z, z, z, z, z, z, z, z, &
      !                  & z, z, z, z, z, z, z, z, z, z, 0.52380952_wp,  c, z, z, z, z, z, z, z, z, &
      !                  & z, z, z, z, z, z, z, z, z, z, -c, 0.52380952_wp, z, z, z, z, z, z, z, z, &
      !                  & z, z, z, z, z, z, z, z, z, z, z, z, 0.61904762_wp,  c, z, z, z, z, z, z, &
      !                  & z, z, z, z, z, z, z, z, z, z, z, z, -c, 0.61904762_wp, z, z, z, z, z, z, &
      !                  & z, z, z, z, z, z, z, z, z, z, z, z, z, z, 0.71428571_wp,  c, z, z, z, z, &
      !                  & z, z, z, z, z, z, z, z, z, z, z, z, z, z, -c, 0.71428571_wp, z, z, z, z, &
      !                  & z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, 0.80952381_wp,  c, z, z, &
      !                  & z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, -c, 0.80952381_wp, z, z, &
      !                  & z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, 0.90476190_wp,  c, &
      !                  & z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, -c, 0.90476190_wp /), &
      !                  & (/ test_size, test_size /), order = (/ 2,1 /))

      allocate(Q(1)) ; !call random_number(Q(1)%data)
      !Q(1)%data = (/ 0.37964151_wp, 1.29390309_wp, 0.34515352_wp, 1.17636033_wp, 0.31379855_wp, 1.06949558_wp, &
      !          & 0.28529197_wp, 0.97233879_wp, 0.25937504_wp, 0.88400808_wp, 0.23581249_wp, 0.80370164_wp, &
      !          & 0.21439044_wp, 0.73069052_wp, 0.19491445_wp, 0.66431200_wp, 0.17720773_wp, 0.60396353_wp, &
      !          & 0.16110955_wp, 0.54909734_wp /)
      Qmat = Q(1)%data
      allocate(Xref(1)) ; call mat_zero(Xref)
      allocate(Xkryl(1:nk)) ; call mat_zero(Xkryl)
      allocate(Xkrylc(1:nk)) ; call mat_zero(Xkrylc)

      Amat = 0.0_wp;
      !forall (i=1:kdim) Amat(i, i) = 1.0_wp
      Amat = A%data
      !do i = 1, kdim
      !   do j = 1, kdim
      !      call random_number(A(i,j))
      !   end do 
      !end do
      tau = 1.0_wp

      !write(*,*) 'Amat'
      !do i = 1, test_size
      !   write(*,'(20F8.4)') Amat(i,1:test_size)
      !enddo
      !write(*,*) 'Qmat test'
      !do i = 1, kdim
      !   write(*,'(1E15.8)') Qmat(i)
      !enddo

      Emat = 0.0_wp
      Xrefmat = 0.0_wp
      call expm(Emat, tau*Amat)
      Xrefmat = matmul(Emat,Qmat)

      !write(*,*) 'expm(Amat)'
      !write(*,*) Xrefmat(1:test_size)
      
      !write(*,*) 'expm(Amat) @ B'
      !do i = 1, test_size
      !   write(*,*) Xrefmat(i), Qmat(i)
      !end do
      Xref(1)%data = Xrefmat

      !write(*,*) 'Emat test'
      !do i = 1, kdim
      !   write(*,'(10E15.8)') Emat(i,1:kdim)
      !enddo

      Emat = 0.0_wp
      do i = 1, nk
         call kexpm(Xkryl(i:i), A, Q(1:1), tau, 1.0e-6_wp, info)
         call Xkryl(i)%axpby(1.0_wp, Xref(1), -1.0_wp)
         difference(i,1) = Xkryl(i)%norm()
         write(*,*) 'norm difference', difference(i,1:1)
         write(*,*) 'y, step',i
         !write(*,*) Xkryl(i)%data
         !write(*,*) Xkrylc(i)%data
      end do


      do i = 1,nk
         write(*,'(a4,"     ",i2,"     ",2E15.8)') 'm = ',i, difference(i,1)
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