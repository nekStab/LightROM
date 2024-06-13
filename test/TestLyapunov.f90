module TestLyapunov
   use LightKrylov
   use LightKrylov, only : wp => dp
   use TestVector
   use TestMatrices
   Use LightROM_LyapunovUtils
   use stdlib_math, only : linspace
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

      ! --> Mesh related parameters.
      real(kind=wp), parameter :: L  = 20.0_wp !> Domain length
      integer      , parameter :: nx = 12      !> Number of grid points (excluding boundaries).
      real(kind=wp), parameter :: dx = L/nx     !> Grid size.

      ! --> Physical parameters.
      complex(kind=wp), parameter :: nu    = cmplx(2.0_wp, 0.2_wp, kind=wp)
      complex(kind=wp), parameter :: gamma = cmplx(1.0_wp, -1.0_wp, kind=wp)
      real(kind=wp)   , parameter :: mu_0  = 0.38_wp
      real(kind=wp)   , parameter :: c_mu  = 0.2_wp
      real(kind=wp)   , parameter :: mu_2  = -0.01_wp
      real(kind=wp)               :: mu(1:nx)

      real(kind=wp)               :: x(1:nx+2)
      
      x = linspace(-L/2, L/2, nx+2)

      x = exp(-(x - 3.0_wp)**2/2.0_wp)

      write(*,*) x


      call check(error, 0.0_wp < rtol)
      
      return
   end subroutine playground

end module TestLyapunov