module LightROM_TestLyapunov
   ! standard library
   use stdlib_math, only : linspace, all_close
   use stdlib_stats_distribution_normal, only: normal => rvs_normal
   use stdlib_linalg, only : svdvals, mnorm
   use stdlib_io_npy, only : save_npy, load_npy
   use stdlib_strings, only : padr
   ! testing library
   use testdrive  , only : new_unittest, unittest_type, error_type, check
   ! LightKrylov for Linear Algebra
   use LightKrylov, only : dp
   use LightKrylov_Logger
   ! LightROM
   use LightROM_Utils
   ! Specific types for testing
   use LightROM_TestUtils
   ! Tests
   Use LightROM_LyapunovUtils  
   use TestUtils
   ! LightControl for reference solution
   use LightControl, only: lyap, dlyap
   
   implicit none
 
   character(len=*), parameter, private :: this_module = 'LightROM_TestUtils'
 
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
            new_unittest("Data POD", test_Proper_Orthogonal_Decomposition_Data_rdp), &
            new_unittest("Balancing Transformation", test_Balancing_Transformation_rdp) &
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
      real(dp) :: norm_direct, norm_LR
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
      norm_direct = mnorm(matmul(Udata, matmul(S, transpose(Udata))) - matmul(Vdata, matmul(G, transpose(Vdata))))
      
      ! Project onto common basis.
      allocate(UTV(ku, kv), VpTV(kv, kv))
      call project_onto_common_basis_rdp(UTV, VpTV, U, V)

      ! Compute frobenius norm from the projection.
      allocate(DLR(ku+kv, ku+kv))
      DLR(    :ku    ,     :ku   ) = S - matmul(UTV,  matmul(G, transpose(UTV)) )
      DLR(ku+1:ku+kv ,     :ku   ) =   - matmul(VpTV, matmul(G, transpose(UTV)) )
      DLR(    :ku    , ku+1:ku+kv) =   - matmul(UTV,  matmul(G, transpose(VpTV)))
      DLR(ku+1:ku+kv , ku+1:ku+kv) =   - matmul(VpTV, matmul(G, transpose(VpTV)))
      norm_LR = mnorm(DLR)

      ! Check correctness.
      err = abs(norm_direct - norm_LR) / max(norm_direct, tiny(1.0_dp))
      call get_err_str(msg, "max err: ", err)
      call check(error, err < rtol_dp)
      call check_test(error, 'test_project_onto_common_basis_rdp', 'Projection consistency', '||X-Y|| = ||X-Y||_LR', msg)
      
      return
   end subroutine test_project_onto_common_basis_rdp

   subroutine test_Proper_Orthogonal_Decomposition_Impulse_rdp(error)
      implicit none
      ! Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      type(state_vector), allocatable :: X0(:)
      class(abstract_vector_rdp), allocatable :: svecs(:)
      type(GL_exponential_prop), allocatable :: prop
      real(dp), dimension(:), allocatable :: svals, sref
      real(dp), dimension(:,:), allocatable :: Q, A

      ! Define test parameters
      real(dp), parameter :: tau = 1.0_dp
      ! Time difference between snapshots
      real(dp), parameter :: Tend = 150.0_dp
      ! Total integration time
      integer :: nprint, i, j, k, ie, is
      integer :: nrank, nstep, nsnap, mode
      real(dp) :: err
      character(len=256) :: msg, info

      integer,  parameter :: irow         = 8
      logical,  parameter :: trans        = .false.
      real(dp), parameter :: tol          = 1e-6_dp
      logical,  parameter :: rescale      = .true.
      integer,  parameter :: rescale_mode = 2
      logical,  parameter :: verbose      = .false.
      
      ! Initialize problem
      call initialize_GL_parameters(X0, A, Q)      
      sref = svdvals(lyap(A, Q))
      if (verbose) then
         nprint = min(8, size(sref))
         do i = 1, ceiling(nprint*1.0_dp/irow)
            is = (i-1)*irow+1; ie = i*irow
            print '(A22,1X,I2,"-",I2,*(1X,F12.8))', padr(' SVD(Xref)',22), is, ie, ( sref(j), j = is, ie )
         end do
         print *, ''
      end if
      
      ! Initialize propagator
      prop = GL_exponential_prop(tau)

      do mode = 1, 2
         ! Compute POD using propagator directly
         call Proper_Orthogonal_Decomposition(svals, prop, X0, tau, Tend, trans, tol, rescale, rescale_mode, svecs=svecs)
         nprint = min(8, size(svals))
         svals(:nprint) = abs(svals(:nprint) - sref(:nprint))
       
         if (verbose) then
               print '(A,I0,A)', 'POD of impulse response, time integration mode ',mode,': Absolute errors in the leading singular values:'
               do i = 1, ceiling(nprint*1.0_dp/irow)
                  is = (i-1)*irow+1; ie = min(i*irow, nprint)
                  print '(1X,A,F6.4,A,I2,A,I2,*(1X,E12.5))', 'SVD err [ dt=', tau,' ]', is, '-', ie, ( svals(j), j = is, ie )
               end do
         end if
       
         err = maxval(svals(:2))
         call get_err_str(msg, "max err: ", err)
         call check(error, err < rtol_dp)
         write(info,'(A,I0,A)') 'test_POD_Imp_', mode, '_rdp'
         call check_test(error, info, 'Leading singular values', 's_1/2 = sPOD_1/2', msg)
      end do
   end subroutine test_Proper_Orthogonal_Decomposition_Impulse_rdp

   subroutine test_Proper_Orthogonal_Decomposition_Data_rdp(error)
      implicit none
      ! Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      type(state_vector), allocatable :: X0(:)
      type(GL_exponential_prop), allocatable :: prop
      real(dp), dimension(:), allocatable :: svals, sref
      real(dp), dimension(:,:), allocatable :: Q, A
      class(abstract_vector_rdp), allocatable :: X(:)   ! Snapshot matrix
      class(abstract_vector_rdp), allocatable :: svecs(:)

      ! Define test parameters
      real(dp), parameter :: tau = 1.0_dp
      ! Time difference between snapshots
      real(dp), parameter :: Tend = 150.0_dp
      ! Total integration time
      integer :: nprint, i, j, k, ie, is
      integer :: nrank, nstep, nsnap, mode
      real(dp) :: err
      character(len=256) :: msg, info

      integer,  parameter :: irow         = 8
      integer,  parameter :: nseries      = 2
      logical,  parameter :: trans        = .false.
      real(dp), parameter :: tol          = 1e-6_dp
      integer,  parameter :: rescale_mode = 2
      logical,  parameter :: verbose      = .false.
      
      ! Initialize problem
      call initialize_GL_parameters(X0, A, Q)
      sref = svdvals(lyap(A, Q))
      
      if (verbose) then
         nprint = min(8, size(sref))
         do i = 1, ceiling(nprint*1.0_dp/irow)
            is = (i-1)*irow+1; ie = i*irow
            print '(A22,1X,I2,"-",I2,*(1X,F12.8))', padr(' SVD(Xref)',22), is, ie, ( sref(j), j = is, ie )
         end do
         print *, ''
      end if

      ! Initialize propagator
      prop = GL_exponential_prop(tau)

      call compute_impulse_response(X, X0, prop, Tend, tau, trans, rescale=.true., rescale_mode=rescale_mode)

      do mode = 1, 2
         ! Compute POD using data matrix
         call Proper_Orthogonal_Decomposition(svals, X, tau, nseries, tol, rescale=.false., svecs=svecs)

         nprint = min(8, size(svals))
         svals(:nprint) = abs(svals(:nprint) - sref(:nprint))
         
         if (verbose) then
            print '(A,I0,A)', 'POD of data matrix, time integration mode ',mode,': Absolute errors in the leading singular values:'
            do i = 1, ceiling(nprint*1.0_dp/irow)
               is = (i-1)*irow+1; ie = min(i*irow, nprint)
               print '(1X,A,F6.4,A,I2,A,I2,*(1X,E12.5))', 'SVD err [ dt=', tau,' ]', is, '-', ie, ( svals(j), j = is, ie )
            end do
         end if

         err = maxval(svals(:2))
         call get_err_str(msg, "max err: ", err)
         call check(error, err < rtol_dp)
         write(info,'(A,I0,A)') 'test_POD_Data_', mode, '_rdp'
         call check_test(error, info, 'Leading singular values', 's_1/2 = sPOD_1/2', msg)
      end do
   end subroutine test_Proper_Orthogonal_Decomposition_Data_rdp

   subroutine test_Balancing_Transformation_rdp(error)
      implicit none
      ! Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      type(state_vector), allocatable :: X0(:), Y0(:)
      type(GL_exponential_prop), allocatable :: prop
      real(dp), dimension(:,:), allocatable :: Q, A
      ! Impule response snapshots
      class(abstract_vector_rdp), allocatable :: X(:), Y(:)   ! Snapshot matrices
      ! Balanced basis
      class(abstract_vector_rdp), allocatable :: T_balanced(:), Tinv_balanced(:), Ttmp(:)
      real(dp), allocatable :: S(:)
      real(dp), dimension(:), allocatable :: Wo(:,:), Wc(:, :)
      ! Reduced operators
      real(dp), allocatable :: Ahat(:,:), Bhat(:,:), Chat(:,:), Xhat(:,:), Yhat(:,:)
      
      ! Define test parameters
      real(dp), parameter :: tau = 1.0_dp
      ! Time difference between snapshots
      real(dp), parameter :: Tend = 100.0_dp
      ! Total integration time
      integer :: i, j, k, rk
      integer :: nrank, nstep, nsnap
      integer :: ndir, nadj, nbal
      real(dp) :: err
      character(len=256) :: msg

      ! Initialize propagator
      prop = GL_exponential_prop(tau)

      ! Initialize forward problem
      call initialize_GL_parameters(X0, A, Q)
      Wc = lyap(A, Q)
      call compute_impulse_response(X, X0, prop, Tend, tau, trans=.false., rescale=.true., rescale_mode=2)

      ! Initialize adjoint problem
      call initialize_GL_parameters(Y0, A, Q, adjoint=.true.)
      Wo = lyap(A, Q)
      call compute_impulse_response(Y, Y0, prop, Tend, tau, trans=.false., rescale=.true., rescale_mode=2)

      call Balancing_Transformation(T_balanced, S, Tinv_balanced, X, Y)

      rk = size(S)
      err = max(norm2(innerprod(Tinv_balanced, T_balanced) - eye(rk)), norm2(innerprod(T_balanced, Tinv_balanced) - eye(rk)))
      call get_err_str(msg, "max err: ", err)
      call check(error, err < rtol_dp)
      call check_test(error, 'test_balancing_transformation', 'Transformation consistency', 'T^{-1} @ T = I', msg)

      !err = norm2(matmul(innerprod(T_balanced, Y), innerprod(Y, T_balanced)) - diag(S))
      !call get_err_str(msg, "max err: ", err)
      !call check(error, err < rtol_dp)
      !call check_test(error, 'test_balancing_transformation', 'Transformation consistency', 'T.T     @ Wo @ T      = Sigma', msg)
      !
      !err = norm2(matmul(innerprod(Tinv_balanced, X), innerprod(X, Tinv_balanced)) - diag(S))
      !call get_err_str(msg, "max err: ", err)
      !call check(error, err < rtol_dp)
      !call check_test(error, 'test_balancing_transformation', 'Transformation consistency', 'T.^{-T} @ Wc @ T^{-1} = Sigma', msg)

      call ROM_Petrov_Galerkin_Projection(Ahat, Bhat, Chat, prop, X0, Y0, T_balanced, Tinv_balanced)
      Xhat = innerprod(Tinv_balanced, X)
      Yhat = innerprod(T_balanced, Y)

      Wc = matmul(Xhat, transpose(Xhat))
      Wo = matmul(Yhat, transpose(Yhat))

      err = norm2(Wc - diag(S))
      call get_err_str(msg, "max err: ", err)
      call check(error, err < rtol_dp)
      call check_test(error, 'test_balancing_transformation', 'Transformation consistency', 'Xhat @ Xhat.T = Sigma', msg)
      err = norm2(Wo - diag(S))
      call get_err_str(msg, "max err: ", err)
      call check(error, err < rtol_dp)
      call check_test(error, 'test_balancing_transformation', 'Transformation consistency', 'Yhat @ Yhat.T = Sigma', msg)

      
   end subroutine test_Balancing_Transformation_rdp

end module LightROM_TestLyapunov
