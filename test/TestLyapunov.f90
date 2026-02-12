module LightROM_TestLyapunov
   ! standard library
   use stdlib_math, only : linspace, all_close
   use stdlib_stats_distribution_normal, only: normal => rvs_normal
   use stdlib_linalg, only : svdvals
   use stdlib_io_npy, only : save_npy, load_npy
   use stdlib_strings, only : padr
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
            new_unittest("project onto common basis", test_project_onto_common_basis_rdp), &
            new_unittest("Impulse Response POD", test_Proper_Orthogonal_Decomposition_Impulse_rdp), &
            new_unittest("Data POD", test_Proper_Orthogonal_Decomposition_Data_rdp) &
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

<<<<<<< HEAD
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
=======
   subroutine test_Proper_Orthogonal_Decomposition_Impulse_rdp(error)
      implicit none
      ! Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      type(state_vector), allocatable :: X0(:)
      type(GL_exponential_prop), allocatable :: prop
      real(dp), dimension(:), allocatable :: svals, sref
      real(dp), dimension(:,:), allocatable :: BBT, A, Xref

      ! Define test parameters
      real(dp), parameter :: tau = 1.0_dp
      ! Time difference between snapshots
      real(dp), parameter :: Tend = 100.0_dp
      ! Total integration time
      integer :: nprint, i, j, k, ie, is
      integer :: nrank, nstep, nsnap
      real(dp) :: res_norm, err
      character(len=256) :: msg

      integer, parameter :: irow = 8
      
      ! Initialize problem
      call initialize_GL_parameters(X0, A, BBT)
      !call solve_lyapunov(Xref, A, BBT)
      call load_npy('test/Xref.npy', Xref)
      res_norm = norm2(matmul(A, Xref) + matmul(Xref, transpose(A)) + BBT)
      !print *, ""
      !print *, 'Residual norm of reference solution: ', res_norm  
      
      sref = svdvals(Xref)
      !nprint = min(8, size(sref))
      !do i = 1, ceiling(nprint*1.0_wp/irow)
      !   is = (i-1)*irow+1; ie = i*irow
      !   print '(A22,1X,I2,"-",I2,*(1X,F12.8))', padr(' SVD(Xref)',22), is, ie, ( sref(j), j = is, ie )
      !end do
      !print *, ''
      
      ! Initialize propagator
      prop = GL_exponential_prop(tau)

      ! Compute POD using propagator directly
      call Proper_Orthogonal_Decomposition(svals, prop, X0, tau, Tend, .false., mode=1)
      nprint = min(8, size(svals))
      svals(:nprint) = (svals(:nprint) - sref(:nprint))**2
      !print *, 'POD of impulse response, time integration mode 1: Absolute errors in the leading singular values:'
      !do i = 1, ceiling(nprint*1.0_wp/irow)
      !   is = (i-1)*irow+1; ie = min(i*irow, nprint)
      !   print '(1X,A,F6.4,A,I2,A,I2,*(1X,E12.5))', 'SVD(XTX) [ dt=', tau,' ]', is, '-', ie, ( svals(j), j = is, ie )
      !end do

      err = maxval(svals(:2))
      call get_err_str(msg, "max err: ", err)
      call check(error, err < rtol_dp)
      call check_test(error, 'test_POD_Imp_1_rdp', info='Leading singular values', eq='s_1/2 = sPOD_1/2', context=msg)

      call Proper_Orthogonal_Decomposition(svals, prop, X0, tau, Tend, .false., mode=2)
      nprint = min(8, size(svals))
      svals(:nprint) = (svals(:nprint) - sref(:nprint))**2
      !print *, 'POD of impulse response, time integration mode 2: Absolute errors in the leading singular values:'
      !do i = 1, ceiling(nprint*1.0_wp/irow)
      !   is = (i-1)*irow+1; ie = min(i*irow, nprint)
      !   print '(1X,A,F6.4,A,I2,A,I2,*(1X,E12.5))', 'SVD(XTX) [ dt=', tau,' ]', is, '-', ie, ( svals(j), j = is, ie )
      !end do

      err = maxval(svals(:2))
      call get_err_str(msg, "max err: ", err)
      call check(error, err < rtol_dp)
      call check_test(error, 'test_POD_Imp_2_rdp', info='Leading singular values', eq='s_1/2 = sPOD_1/2', context=msg)
   end subroutine test_Proper_Orthogonal_Decomposition_Impulse_rdp

   subroutine test_Proper_Orthogonal_Decomposition_Data_rdp(error)
      implicit none
      ! Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      type(state_vector), allocatable :: X0(:)
      type(GL_exponential_prop), allocatable :: prop
      real(dp), dimension(:), allocatable :: svals, sref
      real(dp), dimension(:,:), allocatable :: BBT, A, Xref
      class(state_vector), allocatable :: X(:)   ! Snapshot matrix

      ! Define test parameters
      real(dp), parameter :: tau = 1.0_dp
      ! Time difference between snapshots
      real(dp), parameter :: Tend = 100.0_dp
      ! Total integration time
      integer :: nprint, i, j, k, ie, is
      integer :: nrank, nstep, nsnap
      real(dp) :: res_norm, err
      character(len=256) :: msg

      integer, parameter :: irow = 8
      
      ! Initialize problem
      call initialize_GL_parameters(X0, A, BBT)
      !call solve_lyapunov(Xref, A, BBT)
      call load_npy('test/Xref.npy', Xref)
      res_norm = norm2(matmul(A, Xref) + matmul(Xref, transpose(A)) + BBT)
      !print *, ""
      !print *, 'Residual norm of reference solution: ', res_norm  
      
      sref = svdvals(Xref)
      !nprint = min(8, size(sref))
      !do i = 1, ceiling(nprint*1.0_wp/irow)
      !   is = (i-1)*irow+1; ie = i*irow
      !   print '(A22,1X,I2,"-",I2,*(1X,F12.8))', padr(' SVD(Xref)',22), is, ie, ( sref(j), j = is, ie )
      !end do
      !print *, ''

      ! Initialize propagator
      prop = GL_exponential_prop(tau)

      ! Compute POD using data matrix
      nrank = size(X0)
      nstep = floor(Tend/tau)
      nsnap = nrank*(nstep + 1)
      ! Compute impulse response using propagator
      allocate(X(nsnap))
      k = 1
      do j = 1, nrank ! one series for each initial condition
         call copy(X(k), X0(j))
         do i = 1, nstep ! for the chosen time horizon
            call prop%matvec(X(k), X(k+1))
            k = k + 1
         end do
      end do
      call Proper_Orthogonal_Decomposition(svals, X, tau, nseries=2, mode=1)
      nprint = min(8, size(svals))
      svals(:nprint) = (svals(:nprint) - sref(:nprint))**2
      !print *, 'POD of data matrix, time integration mode 1: Absolute errors in the leading singular values:'
      !do i = 1, ceiling(nprint*1.0_wp/irow)
      !   is = (i-1)*irow+1; ie = min(i*irow, nprint)
      !   print '(1X,A,F6.4,A,I2,A,I2,*(1X,E12.5))', 'SVD(XTX) [ dt=', tau,' ]', is, '-', ie, ( svals(j), j = is, ie )
      !end do

      err = maxval(svals(:2))
      call get_err_str(msg, "max err: ", err)
      call check(error, err < rtol_dp)
      call check_test(error, 'test_POD_Data_1_rdp', info='Leading singular values', eq='s_1/2 = sPOD_1/2', context=msg)

      call Proper_Orthogonal_Decomposition(svals, X, tau, nseries=2, mode=2)
      nprint = min(8, size(svals))
      svals(:nprint) = (svals(:nprint) - sref(:nprint))**2
      !print *, 'POD of data matrix, time integration mode 2: Absolute errors in the leading singular values:'
      !do i = 1, ceiling(nprint*1.0_wp/irow)
      !   is = (i-1)*irow+1; ie = min(i*irow, nprint)
      !   print '(1X,A,F6.4,A,I2,A,I2,*(1X,E12.5))', 'SVD(XTX) [ dt=', tau,' ]', is, '-', ie, ( svals(j), j = is, ie )
      !end do

      err = maxval(svals(:2))
      call get_err_str(msg, "max err: ", err)
      call check(error, err < rtol_dp)
      call check_test(error, 'test_POD_Data_2_rdp', info='Leading singular values', eq='s_1/2 = sPOD_1/2', context=msg)
   end subroutine test_Proper_Orthogonal_Decomposition_Data_rdp

end module LightROM_TestLyapunov
>>>>>>> DLRA
