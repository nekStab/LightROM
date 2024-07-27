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
   use LightKrylov_TestUtils
   use LightKrylov_Utils, only: eigh
   ! LightROM
   use LightROM_Utils
   use LightROM_TestUtils
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
            new_unittest("project onto common basis", test_project_onto_common_basis_rdp), &
            new_unittest("CALE residual norm", test_CALE_res_norm_rdp) &
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
      ! Miscellaneous.
      integer :: kmax, i, j
      real(dp), dimension(:,:), allocatable :: DLR, mu, var
      real(wp) :: norm_direct, norm_LR
      real(dp) :: err, sigma_direct, sigma_projected
      character*256 :: msg

      kmax = max(ku, kv)
      allocate(mu(kmax,kmax), var(kmax, kmax))
      mu = 0.0_dp
      var = 1.0_dp

      ! Initialize bases and coefficients.
      allocate(U(ku), V(kv))
      call init_rand(U); call orthonormalize_basis(U)
      call init_rand(V); call orthonormalize_basis(V)
      allocate(S(ku, ku), G(kv, kv))
      S = normal(mu(:ku,:ku), var(:ku,:ku))
      G = normal(mu(:kv,:kv), var(:kv,:kv))
      
      ! Get data.
      allocate(Udata(test_size, ku), Vdata(test_size, kv))
      call get_data(Udata, U)
      call get_data(Vdata, V)

      ! Compute Frobenius norm directly.
      norm_direct = dense_frobenius_norm(matmul(Udata, matmul(S, transpose(Udata))) - matmul(Vdata, matmul(G, transpose(Vdata))))
      
      ! Project onto common basis.
      allocate(UTV(ku, kv), VpTV(kv, kv))
      call project_onto_common_basis_rdp(UTV, VpTV, U, V)

      ! Compute frobenius norm from the projection.
      allocate(DLR(ku+kv, ku+kv))
      DLR(    :ku    ,     :ku   ) = S - matmul(UTV,  matmul(G, transpose(UTV)) )
      DLR(ku+1:ku+kv ,     :ku   ) =   - matmul(VpTV, matmul(G, transpose(UTV)) )
      DLR(    :ku    , ku+1:ku+kv) =   - matmul(UTV,  matmul(G, transpose(VpTV)))
      DLR(ku+1:ku+kv , ku+1:ku+kv) =   - matmul(VpTV, matmul(G, transpose(VpTV)))
      norm_LR = dense_frobenius_norm(DLR)

      ! Check correctness.
      err = abs(norm_direct - norm_LR)
      call get_err_str(msg, "max err: ", err)
      call check(error, err < rtol_dp)
      call check_test(error, 'test_project_onto_common_basis_rdp', &
                  & info='Equality of difference norm', eq='||X-Y|| = ||X-Y||_LR', context=msg)
      
      return
   end subroutine test_project_onto_common_basis_rdp

   subroutine test_CALE_res_norm_rdp(error)
      ! Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      ! Test linear operator
      type(linop_rdp), allocatable :: A
      ! Test Vectors.
      integer, parameter :: ku = 10
      integer, parameter :: kb = 2
      type(vector_rdp), allocatable :: U(:), B(:)
      ! Coefficient matrices.
      real(dp), allocatable :: S(:, :)
      ! Data matrices.
      real(dp), dimension(:,:), allocatable :: Udata, Bdata, Xdata
      ! Information flag.
      integer :: info
      ! Miscellaneous.
      integer :: i, j
      real(dp), allocatable :: V(:,:), lambda(:)
      real(dp), dimension(test_size, test_size) :: mu, var
      real(dp) :: err, res_direct, res_LR
      type(sym_low_rank_state_rdp) :: X
      character*256 :: msg

      mu = 0.0_dp
      var = 1.0_dp

      ! Initialize basis and coefficients.
      allocate(U(ku), B(kb))
      call init_rand(U); call orthonormalize_basis(U)
      allocate(S(ku, ku), V(ku, ku), lambda(ku))
      S = normal(mu(:ku,:ku), var(:ku,:ku))
      ! ensure symmetry
      S = 0.5*(S + transpose(S))
      ! ensure positive definiteness
      call eigh(S, V, lambda)
      S = matmul(V, matmul(diag(abs(lambda)), transpose(V)))
      ! ensure symmetry which is not exactly assured in the matmul operation
      S = 0.5*(S + transpose(S))
      ! population LR state
      allocate(X%U(ku), source=U)
      allocate(X%S(ku,ku)); X%S = S
      X%rk = ku


      ! Initialize RHS
      call init_rand(B)
     
      ! Initialize operator
      allocate(A)
      A%data = normal(mu, var)

      ! Get data.
      allocate(Udata(test_size, ku))
      call get_data(Udata, U)
      allocate(Bdata(test_size, kb))
      call get_data(Bdata, B)

      ! Compute residual norm directly directly.
      allocate(Xdata(test_size, test_size))
      Xdata = matmul(Udata, matmul(S, transpose(Udata)))
      res_direct = sqrt(sum(svdvals(matmul(A%data, Xdata) + matmul(Xdata, transpose(A%data)) + matmul(Bdata, transpose(Bdata)))**2))
      
      ! Compute low-rank residual norm
      res_LR = CALE_res_norm(X, A, B, 1.0_dp)

      ! Check correctness.
      err = abs(res_LR - res_direct)
      call get_err_str(msg, "max err: ", err)
      call check(error, err < rtol_dp)
      call check_test(error, 'test_CALE_res_norm_rdp', info='Norm equality', eq='||X|| = ||X||_LR', context=msg)
      
      return
   end subroutine test_CALE_res_norm_rdp

end module TestLyapunov
