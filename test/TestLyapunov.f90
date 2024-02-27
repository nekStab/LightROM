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
            new_unittest("Development tests", playground) &
          ]
 
     return
   end subroutine collect_lyapunov_utils_testsuite
   
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
      
      call check(error, 0.0_wp < rtol)
      
      return
   end subroutine playground

end module TestLyapunov