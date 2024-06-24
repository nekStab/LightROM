module TestLyapunov
   ! standard library
   use stdlib_math, only : linspace, all_close
   use stdlib_stats_distribution_normal, only: normal => rvs_normal
   use stdlib_linalg, only : svdvals
   ! testing library
   use testdrive  , only : new_unittest, unittest_type, error_type, check
   ! LightKrylov for Linear Algebra
   use LightKrylov
   use LightKrylov, only : dp, wp => dp
   use LightKrylov_Logger
   use LightKrylov_TestTypes
   ! LightROM
   use LightROM_Utils
   ! Specific types for testing
   use TestUtils
   ! Tests
   Use LightROM_LyapunovUtils
   
   
   implicit none
 
   private :: this_module
   character(len=*), parameter :: this_module = 'LightROM_TestUtils'
 
   public :: collect_lyapunov_utils_testsuite

  contains
 
   !-------------------------------------------
   !-----                                 -----
   !-----     TEST SUITE FOR THE DLRA     -----
   !-----                                 -----
   !-------------------------------------------
 
   subroutine collect_lyapunov_utils_testsuite(testsuite)
     type(unittest_type), allocatable, intent(out) :: testsuite(:)
 
     testsuite = [&
            new_unittest("project onto common basis", test_project_onto_common_basis_rdp) &
          ]
 
     return
   end subroutine collect_lyapunov_utils_testsuite
   
   subroutine test_project_onto_common_basis_rdp(error)
      ! Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      ! Test Vectors.
      integer, parameter :: ku = 10
      integer, parameter :: kv = 15
      type(vector_rdp), allocatable :: U(:), V(:)
      ! Coefficient matrices.
      real(dp), allocatable :: S(:, :), G(:, :)
      ! Data matrices.
      real(dp), allocatable :: Udata(:, :), Vdata(:, :)
      ! Common basis projection results.
      real(dp), allocatable :: UTV(:, :), VpTV(:, :)
      ! Information flag.
      integer :: info
      ! Miscellaneous.
      integer :: kmax
      real(dp) :: mu, var
      real(wp), allocatable :: wrk(:,:)
      real(wp), dimension(:), allocatable :: svals, sdata
      real(dp) :: err, sigma_direct, sigma_projected
      real(dp), dimension(:,:), allocatable :: Ddata, DLR
      character*256 :: msg

      mu = 0.0_dp
      var = 1.0_dp

      ! scratch
      kmax = max(ku,kv)
      allocate(wrk(kmax,kmax));

      ! Initialize bases and coefficients.
      allocate(U(ku), V(kv))
      call init_rand(U); wrk = 0.0_wp; call qr(U, wrk(:ku,:ku), info)
      call init_rand(V); wrk = 0.0_wp; call qr(V, wrk(:kv,:kv), info)
      allocate(S(ku, ku), G(kv, kv))
      S = normal(mu, var); S = 0.5*(S + transpose(S))
      G = normal(mu, var); G = 0.5*(G + transpose(G))
      
      ! Get data.
      allocate(Udata(test_size, ku), Vdata(test_size, kv))
      call get_data(Udata, U)
      call get_data(Vdata, V)

      ! Compute the first singular value directly.
      allocate(Ddata(test_size, test_size)); allocate(sdata(test_size))
      Ddata = matmul(Udata, matmul(S, transpose(Udata))) - matmul(Vdata, matmul(G, transpose(Vdata)))
      sdata = svdvals(Ddata)
      
      ! Project onto common basis.
      allocate(UTV(ku, kv), VpTV(kv, kv))
      call project_onto_common_basis_rdp(UTV, VpTV, U, V)
      call check_info(info, 'project_onto_common_basis_rdp', module=this_module, &
            & procedure='test_project_onto_common_basis_rdp')

      ! Compute the first singular value from the projection.
      allocate(DLR(ku+kv, ku+kv)); allocate(svals(ku+kv))
      DLR(    :ku    ,     :ku   ) = S - matmul(UTV,  matmul(G, transpose(UTV)) )
      DLR(ku+1:ku+kv ,     :ku   ) =   - matmul(VpTV, matmul(G, transpose(UTV)) )
      DLR(    :ku    , ku+1:ku+kv) =   - matmul(UTV,  matmul(G, transpose(VpTV)))
      DLR(ku+1:ku+kv , ku+1:ku+kv) =   - matmul(VpTV, matmul(G, transpose(VpTV)))
      svals = svdvals(DLR)

      ! Check correctness.
      err = abs(sdata(1) - svals(1))
      call get_err_str(msg, "max err: ", err)
      call check(error, err < rtol_dp)
      call check_test(error, 'test_project_onto_common_basis_rdp', &
                  & info='Singular value comparison', eq='s_1 = s(LR)_1', context=msg)
      
      return
   end subroutine test_project_onto_common_basis_rdp

end module TestLyapunov