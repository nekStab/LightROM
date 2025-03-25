module LightROM_TestLyapunov
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
   ! LightROM
   use LightROM_Utils
   ! Specific types for testing
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
            new_unittest("Impulse Response POD", test_Proper_Orthogonal_Decomposition_Impulse_rdp) &
            !new_unittest("project onto common basis", test_project_onto_common_basis_rdp) &
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

   subroutine test_Proper_Orthogonal_Decomposition_Impulse_rdp(error)
      ! Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      type(state_vector), allocatable :: X0(:)
      type(GL_exponential_prop), allocatable :: prop
      real(dp), dimension(:), allocatable :: svals
      real(dp), dimension(:,:), allocatable :: BBT, A, Xref
      class(state_vector), allocatable :: X(:)   ! Snapshot matrix

      ! Define test parameters
      real(dp), parameter :: tau = 1.0_dp
      ! Time difference between snapshots
      real(dp), parameter :: Tend = 100.0_dp
      ! Total integration time
      integer :: nprint, i, j, k, irow, ie, is
      integer :: nrank, nstep, nsnap

      irow = 8

      ! Initialize problem
      call initialize_GL_parameters(X0, A, BBT)
      call solve_lyapunov(Xref, A, BBT)
      
      svals = svdvals(Xref)
      nprint = min(8, size(svals))
      do i = 1, ceiling(nprint*1.0_wp/irow)
         is = (i-1)*irow+1; ie = i*irow
         print '(A22,1X,I2,"-",I2,*(1X,F16.12))', 'SVD(X_BS)            ', is, ie, ( svals(j), j = is, ie )
      end do
      print *, ''
      
      ! Initialize propagator
      prop = GL_exponential_prop(tau)

      ! Compute POD using propagator directly
      call Proper_Orthogonal_Decomposition(svals, prop, X0, tau, Tend, .false., mode=1)
      nprint = min(8, size(svals))
      do i = 1, ceiling(nprint*1.0_wp/irow)
         is = (i-1)*irow+1; ie = min(i*irow, nprint)
         print '(1X,A,F6.4,A,I2,A,I2,*(1X,F16.12))', 'SVD(XTX) [ dt=', tau,' ]', is, '-', ie, ( svals(j), j = is, ie )
      end do
      !call Proper_Orthogonal_Decomposition(svals, prop, X0, tau, Tend, .false., mode=2)
      !do i = 1, ceiling(nprint*1.0_wp/irow)
      !   is = (i-1)*irow+1; ie = min(i*irow, nprint)
      !   print '(1X,A,F6.4,A,I2,A,I2,*(1X,F16.12))', 'SVD(XTX) [ dt=', tau,' ]', is, '-', ie, ( svals(j), j = is, ie )
      !end do
      !call Proper_Orthogonal_Decomposition(svals, prop, X0, tau, Tend, .false., mode=3)
      !nprint = min(8, size(svals))
      !do i = 1, ceiling(nprint*1.0_wp/irow)
      !   is = (i-1)*irow+1; ie = min(i*irow, nprint)
      !   print '(1X,A,F6.4,A,I2,A,I2,*(1X,F16.12))', 'SVD(XTX) [ dt=', tau,' ]', is, '-', ie, ( svals(j), j = is, ie )
      !end do

      ! Compute POD using data matrix
      nrank = size(X0)
      nstep = floor(Tend/tau)
      nsnap = nrank*(nstep + 1)
      allocate(X(nsnap))
      ! Compute impulse response using propagator
      k = 1
      do j = 1, 2 ! for each initial condition
         call copy(X(k), X0(j))
         do i = 1, nstep ! for the chosen time horizon
            call prop%matvec(X(k), X(k+1))
            k = k + 1
         end do
      end do
      call Proper_Orthogonal_Decomposition(svals, X, tau, mode=1)
      nprint = min(8, size(svals))
      do i = 1, ceiling(nprint*1.0_wp/irow)
         is = (i-1)*irow+1; ie = min(i*irow, nprint)
         print '(1X,A,F6.4,A,I2,A,I2,*(1X,F16.12))', 'SVD(XTX) [ dt=', tau,' ]', is, '-', ie, ( svals(j), j = is, ie )
      end do
      !call Proper_Orthogonal_Decomposition(svals, X, tau, mode=2)
      !nprint = min(8, size(svals))
      !do i = 1, ceiling(nprint*1.0_wp/irow)
      !   is = (i-1)*irow+1; ie = min(i*irow, nprint)
      !   print '(1X,A,F6.4,A,I2,A,I2,*(1X,F16.12))', 'SVD(XTX) [ dt=', tau,' ]', is, '-', ie, ( svals(j), j = is, ie )
      !end do
      !call Proper_Orthogonal_Decomposition(svals, X, tau, mode=3)
      !nprint = min(8, size(svals))
      !do i = 1, ceiling(nprint*1.0_wp/irow)
      !   is = (i-1)*irow+1; ie = min(i*irow, nprint)
      !   print '(1X,A,F6.4,A,I2,A,I2,*(1X,F16.12))', 'SVD(XTX) [ dt=', tau,' ]', is, '-', ie, ( svals(j), j = is, ie )
      !end do
   end subroutine test_Proper_Orthogonal_Decomposition_Impulse_rdp

end module LightROM_TestLyapunov