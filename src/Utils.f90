module LightROM_Utils
   ! stdlib
   use stdlib_strings, only: padl
   use stdlib_linalg, only : eye, diag, svd, svdvals, is_symmetric
   use stdlib_optval, only : optval
   use stdlib_stats_distribution_normal, only: normal => rvs_normal
   use stdlib_logger, only : logger => global_logger
   ! LightKrylov for Linear Algebra
   use LightKrylov
   use LightKrylov, only : dp
   use LightKrylov_Logger, only: log_message, log_information, log_warning, log_debug, check_info, stop_error
   use LightKrylov_AbstractVectors
   use LightKrylov_BaseKrylov, only : orthogonalize_against_basis
   use LightKrylov_Utils, only : abstract_opts, sqrtm
   ! LightROM
   use LightROM_AbstractLTIsystems
   
   implicit none 
   character(len=*), parameter, private :: this_module = 'LR_Utils'

   public :: dlra_opts
   public :: coefficient_matrix_norm, increment_norm, low_rank_CALE_residual_norm
   public :: is_converged
   public :: project_onto_common_basis
   public :: Balancing_Transformation
   public :: LQR_gain
   public :: LQE_gain
   public :: ROM_Petrov_Galerkin_Projection
   public :: ROM_Galerkin_Projection
   public :: Proper_Orthogonal_Decomposition

   interface Balancing_Transformation
      module procedure Balancing_Transformation_rdp
   end interface

   interface LQR_gain
      module procedure LQR_gain_vector_rdp
      module procedure LQR_gain_matrix_rdp
   end interface

   interface LQE_gain
      module procedure LQE_gain_vector_rdp
      module procedure LQE_gain_matrix_rdp
   end interface

   interface ROM_Petrov_Galerkin_Projection
      module procedure ROM_Petrov_Galerkin_Projection_rdp
   end interface

   interface ROM_Galerkin_Projection
      module procedure ROM_Galerkin_Projection_rdp
   end interface

   interface Proper_Orthogonal_Decomposition
      module procedure Proper_Orthogonal_Decomposition_Impulse_rdp
      module procedure Proper_Orthogonal_Decomposition_Data_rdp
   end interface

   interface rescale_snapshots
      module procedure rescale_snapshots_rdp
   end interface

   interface project_onto_common_basis
      module procedure project_onto_common_basis_rdp
   end interface

   type, extends(abstract_opts), public :: dlra_opts
      !! Options container for the (rank-adaptive) projector-splitting dynalical low-rank approximation
      !
      ! BASIC OPTIONS
      !
      !! integrator
      integer :: mode = 1
      !! Time integration mode. Only 1st order (Lie splitting - mode 1) and 
      !! 2nd order (Strang splitting - mode 2) are implemented. (default: 1)
      real(dp) :: tau = 1.0_dp
      !! Requested integration step
      real(dp) :: Tend = 1.0_dp
      !! Final integration time
      integer :: nsteps = 0
      !! Total number of timesteps
      !
      ! CONVERGENCE CHECK
      !
      integer :: chkstep = 10
      !! Time step interval at which convergence is checked and runtime information is printed (default: 10)
      real(dp) :: chktime = 1.0_dp
      !! Simulation time interval at which convergence is checked and runtime information is printed (default: 1.0)
      logical :: chkctrl_time = .true.
      !! Use time instead of timestep control (default: .true.)
      real(dp) :: inc_tol = 1e-6_dp
      !! Tolerance on the increment for convergence (default: 1e-6)
      logical :: relative_inc = .true.
      !! Tolerance control: Use relative values for convergence (default = .true.)
      !
      ! RANK-ADPATIVE SPECIFICS
      !
      logical :: if_rank_adaptive = .true.
      !! Allow rank-adaptivity
      integer :: rk_init = 1
      !! Guess for the initial rank
      integer :: n_init = 5
      !! Number of steps to determine initial rank
      real(dp) :: tol = 1e-6_dp
      !! Tolerance on the extra singular value to determine rank-adaptation
      integer :: err_est_step = 10
      !! Time step interval for recomputing the splitting error estimate (only of use_err_est = .true.)
      integer :: rk_reduction_lock = 0
      !! Current value of the rank reduction lock
      integer :: rk_reduction_barrier = 10
      !! Reset value of the rank reduction lock
   contains
      procedure, pass(self), public :: init => dlra_opts_initialize
   end type

contains

   subroutine dlra_opts_initialize(self)
      class(dlra_opts), intent(inout) :: self
      ! internal
      character(len=*), parameter :: this_procedure = 'dlra_opts_initialize'
      character(len=256) :: msg
      type(dlra_opts) :: opts_default
      integer :: itmp
      real(dp) :: tau_eff
      real(dp), parameter :: tol = 1e-7

      opts_default = dlra_opts()
      !
      ! INTEGRATION TIME
      !
      itmp = nint(self%Tend/self%tau)
      if (self%Tend < 0.0_dp .or. self%tau < 0.0_dp .or. itmp < 1) then
         write(msg,'(2(A,E12.5))') 'Invalid Tend/tau combination specified. Tend = ', self%Tend, ', tau = ', self%tau
         call stop_error(msg, this_module, this_procedure)
      end if 
      if (self%nsteps /= 0 .and. self%nsteps /= itmp) then
         write(msg,'(2(A,E12.5))') 'Specified nsteps does not match Tend/tau: Tend = ', self%Tend, ', tau = ', self%tau
         call log_warning(msg, this_module, this_procedure)
         write(msg,'(A,I0)') 'Resetting. nsteps = ', itmp
         call log_warning(msg, this_module, this_procedure)
      end if
      self%nsteps = itmp
      tau_eff = self%Tend/self%nsteps
      if (abs(self%tau-tau_eff) < tol) then
         write(msg,'(3(A,E12.5))') 'Timestep reset to match specified Tend = ', self%Tend, ': tau = ', self%tau, ' -> ', tau_eff
         call log_warning(msg, this_module, this_procedure)
         self%tau = tau_eff
      end if
      !
      ! TEMPORAL ORDER
      !
      if ( self%mode > 2 ) then
         write(msg,'(A)') "Time-integration order for the operator splitting of d > 2 &
                      & requires adjoint solves and is not implemented. Resetting torder = 2." 
         call log_message(msg, this_module, this_procedure)
      else if ( self%mode < 1 ) then
         write(msg,'(A,I0)') "Invalid time-integration order specified: ", self%mode
         call stop_error(msg, this_module, this_procedure)
      endif
      !
      ! CONVERGENCE CHECK
      !
      if (self%chkctrl_time) then
         if (self%chktime <= 0.0_dp) then
            self%chktime = opts_default%chktime
            write(msg,'(A,E12.5,A)') 'Invalid chktime. Reset to default (',  self%chktime,')'
            call log_warning(msg, this_module, this_procedure)
         end if
         self%chkstep = max(1, NINT(self%chktime/self%tau))
         write(msg,'(A,E12.5,A,I0,A)') 'Convergence check every ', self%chktime, ' time units (', self%chkstep, ' steps)'
         call log_information(msg, this_module, this_procedure)
      else
         if (self%chkstep <= 0) then
            self%chkstep = opts_default%chkstep
            write(msg,'(A,I0,A)') "Invalid chktime. Reset to default (",  self%chkstep,")"
            call log_warning(msg, this_module, this_procedure)
         end if
         self%chktime = self%tau*self%chkstep
         write(msg,'(A,I0,A)') 'Convergence check every ', self%chkstep, ' steps (based on steps).'
         call log_information(msg, this_module, this_procedure)
      end if
      
   end subroutine dlra_opts_initialize

   subroutine Balancing_Transformation_rdp(T, S, Tinv, Xc, Yo)
      !! Computes the the biorthogonal balancing transformation \( \mathbf{T}, \mathbf{T}^{-1} \) from the
      !! low-rank representation of the SVD of the controllability and observability Gramians, \( \mathbf{W}_c \) 
      !! and \( \mathbf{W}_o \) respectively, given as:
      !! \[ \begin{align}
      !!    \mathbf{W}_c &= \mathbf{X}_c \mathbf{X}_c^T \\
      !!    \mathbf{W}_o &= \mathbf{Y}_o \mathbf{Y}_o^T
      !! \end{align} \]
      !!
      !! Given the SVD of the cross-Gramian:
      !! $$ \mathbf{X}_c^T \mathbf{Y}_o = \mathbf{U} \mathbf{S} \mathbf{V}^T $$
      !! the balancing transformation and its inverse are given by:
      !! \[ \begin{align}
      !!            \mathbf{T}   &= \mathbf{X}_o \mathbf{S}_o^{1/2} \mathbf{V} \mathbf{S}^{-1/2} \\
      !!            \mathbf{Tinv}^T &= \mathbf{Y}_c \mathbf{S}_c^{1/2} \mathbf{U} \mathbf{S}^{-1/2} 
      !! \end{align} \]
      !! Note: In the current implementation, the numerical rank of the SVD is not considered.
      class(abstract_vector_rdp),          intent(out)   :: T(:)
      !! Balancing transformation
      real(dp),                            intent(out)   :: S(:)
      !! Singular values of the BT
      class(abstract_vector_rdp),          intent(out)   :: Tinv(:)
      !! Inverse balancing transformation
      class(abstract_vector_rdp),          intent(in)    :: Xc(:)
      !! Low-rank representation of the Controllability Gramian
      class(abstract_vector_rdp),          intent(in)    :: Yo(:)
      !! Low-rank representation of the Observability Gramian

      ! internal variables
      integer                                :: i, rkc, rko, rk, rkmin
      real(dp),                  allocatable :: LRCrossGramian(:,:)
      real(dp),                  allocatable :: Swrk(:,:)
      real(dp),                  allocatable :: Sigma(:)
      real(dp),                  allocatable :: V(:,:), W(:,:)

      rkc   = size(Xc)
      rko   = size(Yo)
      rk    = max(rkc, rko)
      rkmin = min(rkc, rko) 

      ! compute inner product with Gramian bases and compte SVUD
      allocate(LRCrossGramian(rkc,rko)); allocate(V(rko,rko)); allocate(W(rkc,rkc))
      LRCrossGramian = innerprod(Xc, Yo)
      call svd(LRCrossGramian, S, V, W)

      allocate(Sigma(rkmin))
      do i = 1, rkmin
         Sigma(i) = 1/sqrt(S(i))
      enddo
      block
         class(abstract_vector_rdp), allocatable :: Xwrk(:)
         call linear_combination(Xwrk, Yo(1:rkmin), matmul(W(1:rkmin,1:rkmin), diag(Sigma)))
         call copy(T(1:rkmin), Xwrk)
         call linear_combination(Xwrk, Xc(1:rkmin), matmul(V(1:rkmin,1:rkmin), diag(Sigma)))
         call copy(Tinv(1:rkmin), Xwrk)
      end block
         
   end subroutine Balancing_Transformation_rdp

   subroutine LQR_gain_vector_rdp(KT, X, B, Rinv)
      class(abstract_vector_rdp), allocatable, intent(out) :: KT
      !! LGR gains
      class(abstract_sym_low_rank_state_rdp), intent(in) :: X
      !! Low rank solution of current solution
      class(abstract_vector_rdp), intent(in) :: B
      !! System inputs
      real(dp), intent(in) :: Rinv
      !! Inverse control cost

      ! internal variables
      real(dp), allocatable :: proj(:), wrk(:)

      ! K = R^{-1} @ B.T @ Pbut we compute
      ! K.T = P @ B @ R^{-1} = U @ S @ U.T @ B @ R^{-1}
      proj = innerprod(X%U(:X%rk), B)
      wrk  = matmul(X%S(:X%rk,:X%rk), proj * Rinv)
      call linear_combination(KT, X%U(:X%rk), wrk)
   end subroutine LQR_gain_vector_rdp

   subroutine LQR_gain_matrix_rdp(KT, X, B, Rinv)
      class(abstract_vector_rdp), allocatable, intent(out) :: KT(:)
      !! LGR gains
      class(abstract_sym_low_rank_state_rdp), intent(in) :: X
      !! Low rank solution of current solution
      class(abstract_vector_rdp), intent(in) :: B(:)
      !! System inputs
      real(dp), intent(in) :: Rinv(:,:)
      !! Inverse control cost

      ! internal variables
      real(dp), allocatable :: proj(:,:), wrk(:,:)

      ! K = R^{-1} @ B.T @ P but we compute
      ! K.T = P @ B @ R^{-1}
      proj = innerprod(X%U(:X%rk), B)
      wrk  = matmul(X%S(:X%rk,:X%rk), matmul(proj, Rinv))
      call linear_combination(KT, X%U(:X%rk), wrk)  
      
   end subroutine LQR_gain_matrix_rdp
   
   subroutine LQE_gain_vector_rdp(L, X, CT, Vinv)
      class(abstract_vector_rdp), allocatable, intent(out) :: L
      !! Kalman gains
      class(abstract_sym_low_rank_state_rdp), intent(in) :: X
      !! Low rank solution of current solution
      class(abstract_vector_rdp), intent(in) :: CT
      !! Sensors
      real(dp), intent(in) :: Vinv
      !! Inverse sensor noise variance

      ! internal variables
      real(dp), allocatable :: proj(:)

      ! L = P @ C.T @ V^{-1}
      allocate(L, source=CT); call L%scal(Vinv)
      proj = innerprod(X%U(:X%rk), L)
      call linear_combination(L, X%U(:X%rk), matmul(X%S(:X%rk,:X%rk), proj))
      
   end subroutine LQE_gain_vector_rdp
   
   subroutine LQE_gain_matrix_rdp(L, X, CT, Vinv)
      class(abstract_vector_rdp), allocatable, intent(out) :: L(:)
      !! Kalman gains
      class(abstract_sym_low_rank_state_rdp), intent(in) :: X
      !! Low rank solution of current solution
      class(abstract_vector_rdp), intent(in) :: CT(:)
      !! Sensors
      real(dp), intent(in) :: Vinv(:,:)
      !! Inverse sensor noise variance

      ! internal variables
      real(dp), allocatable :: proj(:,:)

      ! L = P @ C.T @ V^{-1}
      call linear_combination(L, CT, Vinv)
      proj = innerprod(X%U(:X%rk), L)
      call linear_combination(L, X%U(:X%rk), matmul(X%S(:X%rk,:X%rk), proj))
      
   end subroutine LQE_gain_matrix_rdp

   subroutine ROM_Petrov_Galerkin_Projection_rdp(Ahat, Bhat, Chat, D, LTI, T, Tinv)
      !! Computes the Reduced-Order Model of the input LTI dynamical system via Petrov-Galerkin projection 
      !! using the biorthogonal projection bases \( \mathbf{V} \) and \( \mathbf{W} \) with 
      !! \( \mathbf{W}^T \mathbf{V} = \mathbf{I} \).
      !! 
      !! Given an LTI system defined by the matrices \( \mathbf{A}, \mathbf{B}, \mathbf{C}, \mathbf{D}\), 
      !! the matrices \( \hat{\mathbf{A}}, \hat{\mathbf{B}}, \hat{\mathbf{C}}, \hat{\mathbf{D}}\) of the 
      !! projected LTI system are given by:
      !! \[
      !!     \hat{\mathbf{A}} = \mathbf{W}^T \mathbf{A} \mathbf{V}, \qquad
      !!     \hat{\mathbf{B}} = \mathbf{W}^T \mathbf{B}, \qquad
      !!     \hat{\mathbf{C}} = \mathbf{C} \mathbf{V}, \qquad
      !!     \hat{\mathbf{D}} = \mathbf{D} .
      !! \]
      real(dp),            allocatable, intent(out)    :: Ahat(:, :)
      !! Reduced-order dynamics matrix.
      real(dp),            allocatable, intent(out)    :: Bhat(:, :)
      !! Reduced-order input-to-state matrix.
      real(dp),            allocatable, intent(out)    :: Chat(:, :)
      !! Reduced-order state-to-output matrix.
      real(dp),            allocatable, intent(out)    :: D(:, :)
      !! Feed-through matrix
      class(abstract_lti_system_rdp),   intent(inout)  :: LTI
      !! Large-scale LTI to project
      class(abstract_vector_rdp),       intent(in)     :: T(:)
      !! Balancing transformation
      class(abstract_vector_rdp),       intent(in)     :: Tinv(:)
      !! Inverse balancing transformation

      ! internal variables
      integer                                          :: i, rk, rkc, rkb
      class(abstract_vector_rdp),       allocatable    :: Uwrk(:)
      real(dp),                         allocatable    :: Cwrk(:, :)

      rk  = size(T)
      rkb = size(LTI%B)
      rkc = size(LTI%CT)
      allocate(Uwrk(rk), source=T(1)); call zero_basis(Uwrk)
      allocate(Ahat(1:rk, 1:rk ));                  Ahat = 0.0_dp
      allocate(Bhat(1:rk, 1:rkb));                  Bhat = 0.0_dp
      allocate(Cwrk(1:rk, 1:rkc));                  Cwrk = 0.0_dp
      allocate(Chat(1:rkc,1:rk ));                  Chat = 0.0_dp
      allocate(D(1:size(LTI%D,1),1:size(LTI%D,2))); D    = 0.0_dp

      do i = 1, rk
         call LTI%A%matvec(Tinv(i), Uwrk(i))
      end do
      Ahat = innerprod(T, Uwrk)
      Bhat = innerprod(T, LTI%B)
      Cwrk = innerprod(LTI%CT, Tinv)
      Chat = transpose(Cwrk)
      D = LTI%D

   end subroutine ROM_Petrov_Galerkin_Projection_rdp

   subroutine ROM_Galerkin_Projection_rdp(Ahat, Bhat, Chat, D, LTI, T)
      !! Computes the Reduced-Order Model of the input LTI dynamical system via Galerkin projection using 
      !! the orthogonal projection basis \( \mathbf{V} \) with \( \mathbf{V}^T \mathbf{V} = \mathbf{I} \).
      !! 
      !! Given an LTI system defined by the matrices \( \mathbf{A}, \mathbf{B}, \mathbf{C}, \mathbf{D}\), 
      !! the matrices \( \hat{\mathbf{A}}, \hat{\mathbf{B}}, \hat{\mathbf{C}}, \hat{\mathbf{D}}\) of the projected LTI system is given by:
      !! \[
      !!     \hat{\mathbf{A}} = \mathbf{V}^T \mathbf{A} \mathbf{V}, \qquad
      !!     \hat{\mathbf{B}} = \mathbf{V}^T \mathbf{B}, \qquad
      !!     \hat{\mathbf{C}} = \mathbf{C} \mathbf{V}, \qquad
      !!     \hat{\mathbf{D}} = \mathbf{D} .
      !! \]
      real(dp),            allocatable, intent(out)    :: Ahat(:, :)
      !! Reduced-order dynamics matrix.
      real(dp),            allocatable, intent(out)    :: Bhat(:, :)
      !! Reduced-order input-to-state matrix.
      real(dp),            allocatable, intent(out)    :: Chat(:, :)
      !! Reduced-order state-to-output matrix.
      real(dp),            allocatable, intent(out)    :: D(:, :)
      !! Feed-through matrix
      class(abstract_lti_system_rdp),   intent(inout)  :: LTI
      !! Large-scale LTI to project
      class(abstract_vector_rdp),       intent(inout)  :: T(:)
      !! Balancing transformation

      call ROM_Petrov_Galerkin_Projection(Ahat, Bhat, Chat, D, LTI, T, T)

   end subroutine ROM_Galerkin_Projection_rdp

   subroutine Proper_Orthogonal_Decomposition_Impulse_rdp(svals, prop, X0, tau, Tend, trans, mode, svecs)
      !! Computes the Proper Orthogonal Decomposition (POD) of the impulse response to the input vector based on the
      !! exponential propagator prop.
      !! 
      !! The POD is computed using the singular value decomposition of the weighted inner product of the snapshot matrix
      !! \[
      !!     (\lambda_i, \xi_i) = \mbox{svd}(X^T W X)
      !! \]
      !! where the snapshot matrix \( X \) is given by
      !! \[
      !!    X = [ \: x\sqrt{\tau} \: ]
      !! \]
      !! and the individual snapshots \( x_i \) are recursively computed as \( x_{i+1} = e^{\tau A} x_i \) using the 
      !! exponential propagator for each initial condition. The integration weights \( W \) are automatically applied 
      !! in the computation of the inner product.
      !!
      !! NOTE: 
      !!    1. We use the snapshot method which assumes that the number of snapshots is smaller than the number of d
      !!       egrees of freedom of the problem
      !!    2. We store all the snapshots in memory before performing the inner product.
      !!    3. We assume that the integration time tau is already set for the propagator
      !!
      real(dp), allocatable, intent(out) :: svals(:)
      !! POD singular values
      class(abstract_linop_rdp), intent(inout) :: prop
      !! Exponential propagator
      class(abstract_vector_rdp), intent(in) :: X0(:)
      !! Initial condition (impulse)
      real(dp), intent(in) :: tau 
      !! Integration time between subsequent snapshots
      real(dp), intent(in) :: Tend
      !! Time horizon for the POD computation
      logical, optional, intent(in) :: trans
      !! Direct of adjoint mode (default: direct)
      integer, optional, intent(in) :: mode
      !! Time integration mode
      class(abstract_vector_rdp), optional, allocatable, intent(out) :: svecs(:)
      !! Singular vectors

      ! internals
      class(abstract_vector_rdp), allocatable :: X(:)   ! Snapshot matrix
      real(dp), allocatable :: XTX(:,:)  ! Inner product matrix
      real(dp), allocatable :: U(:,:), VT(:,:) ! singular vectors
      integer :: i, j, k
      integer :: nsnap, nstep, nrank
      logical :: transpose

      transpose = optval(trans, .false.)

      ! Data sizes
      nrank = size(X0)
      nstep = floor(Tend/tau)
      nsnap = nrank*(nstep + 1)
      
      ! Compute impulse response
      allocate(X(nsnap), source=X0(1)) ; call zero_basis(X)
      k = 1
      do j = 1, nrank ! one series for each initial condition
         call copy(X(k), X0(j))
         do i = 1, nstep ! for the chosen time horizon
            if (transpose) then
               call prop%rmatvec(X(k), X(k+1))
            else
               call prop%matvec(X(k), X(k+1))
            end if
            k = k + 1
         end do
      end do

      ! Rescale response for POD
      do j = 1, nrank
         call rescale_snapshots(X((j-1)*nstep+1:j*nstep), tau, mode)
      end do
      
      ! Compute cross-correlation
      allocate(XTX(nsnap,nsnap))
      XTX = gram(X)
      if (.not. present(svecs)) then
         ! Compute only the POD singular values
         svals = svdvals(XTX)
      else
         ! Compute POD singular values and vectors
         call svd(XTX, svals, U, VT)
         ! Project data matrix onto principal axes
         call linear_combination(svecs, X, U)
      end if
   end subroutine Proper_Orthogonal_Decomposition_Impulse_rdp

   subroutine Proper_Orthogonal_Decomposition_Data_rdp(svals, X, tau, nseries, mode, svecs)
      !! Computes the Proper Orthogonal Decomposition (POD) of the input vector of data snapshots taken at constant time
      !! intervals tau
      !! 
      !! The POD is computed using the singular value decomposition of the weighted inner product of the snapshot matrix
      !! \[
      !!     (\lambda_i, \xi_i) = \mbox{svd}(X^T W X)
      !! \]
      !! where the snapshot matrix \( X \) is given by
      !! \[
      !!    X = [ \: x\sqrt{\tau} \: ]
      !! \]
      !! The integration weights \( W \) are automatically applied in the computation of the inner product.
      !!
      !! NOTE: 
      !!    1. We use the snapshot method which assumes that the number of snapshots is smaller than the number of
      !!       degrees of freedom of the problem
      !!    2. The input data matrix is destroyed.
      !!    3. We assume that the the snapshots are equidistant in time
      !!
      real(dp), allocatable, intent(out) :: svals(:)
      !! POD singular values
      class(abstract_vector_rdp), intent(inout) :: X(:)
      !! Snapshots
      real(dp), intent(in) :: tau 
      !! Integration time between subsequent snapshots
      integer, optional, intent(in) :: nseries
      !! Number of independent time series in the data
      integer, optional, intent(in) :: mode
      !! Time integration mode
      class(abstract_vector_rdp), optional, allocatable, intent(out) :: svecs(:)
      !! Singular vectors

      ! internals
      real(dp), allocatable :: XTX(:,:)  ! Inner product matrix
      real(dp), allocatable :: U(:,:), VT(:,:) ! singular vectors
      integer :: j
      integer :: nsnap, nstep, nrank
      character(len=128) :: msg

      ! Data sizes
      nrank = optval(nseries, 1)
      nsnap = size(X)
      nstep = nsnap/nrank - 1
      if (nsnap .ne. (nstep+1)*nrank) then
         write(msg,'(2(A,I0,A))') "Input data of size ", nsnap, " cannot be partitioned into ", nrank, " series."
         call stop_error(msg, this_module, 'Proper_Orthogonal_Decomposition_Data_rdp')
      end if

      ! Rescale response for POD
      do j = 1, nrank
         call rescale_snapshots(X((j-1)*nstep+1:j*nstep), tau, mode)
      end do
      
      ! Compute cross-correlation
      allocate(XTX(nsnap,nsnap))
      XTX = gram(X)
      ! Compute POD
      if (.not. present(svecs)) then
         ! Compute only the POD singular values
         svals = svdvals(XTX)
      else
         ! Compute POD singular values and vectors
         call svd(XTX, svals, U, VT)
         ! Project data matrix onto principal axes
         call linear_combination(svecs, X, U)
      end if
   end subroutine Proper_Orthogonal_Decomposition_Data_rdp

   subroutine rescale_snapshots_rdp(X, tau, mode)
      !! Apply integration weights to the columns of a data matrix for temporal integration. 
      !! The columns are assumed to be equispaced snapshots at intervals `tau` based on the integration method `mode`.
      class(abstract_vector_rdp), intent(inout) :: X(:)
      !! Snapshot matrix to be rescaled
      real(dp), intent(in) :: tau
      !! Time difference between snapshots in X (assumed constant)
      integer, optional, intent(in) :: mode
      !! Temporal integration mode
      ! internal
      integer :: i, n, mode_
      real(dp) :: w
      character(len=128) :: msg
      mode_ = optval(mode, 1)
      n = size(X)
      ! Sanity checks
      if (tau <= 0.0_dp) then
         write(msg,'(A,F16.12)') "Invalid integration step tau= ", tau
         call stop_error(msg, this_module, 'rescale_snapshots')
      end if
      select case (mode_)
      case (1)
         ! uniform weights
         w = tau
         do i = 1, n
            call X(i)%scal(sqrt(w))
         end do
      case (2)
         ! Trapezoid rule
         w = tau
         call X(1)%scal(sqrt(0.5_dp*w))
         do i = 2, n - 1
            call X(i)%scal(sqrt(w))
         end do
         call X(n)%scal(sqrt(0.5_dp*w))
      case default
         write(msg,'(A,I0)') "The selected integration method is not implemented: ", mode_
         call stop_error(msg, this_module, 'rescale_snapshots')
      end select
   end subroutine rescale_snapshots_rdp

   subroutine project_onto_common_basis_rdp(UTV, VpTV, U, V)
      !! Computes the common orthonormal basis of the space spanned by the union of the input Krylov bases 
      !! \( [ \mathbf{U}, \mathbf{V} ] \) by computing \( \mathbf{V_\perp} \) as an orthonormal basis of 
      !! \( \mathbf{V} \) lying in the orthogonal complement of \( \mathbf{U} \) given by
      !! \[
      !!    \mathbf{V_\perp}, R = \text{qr}( \mathbf{V} - \mathbf{U} \mathbf{U}^T \mathbf{V} )
      !! \[
      !!
      !! NOTE: The orthonormality of \( \mathbf{U} \) is assumed and not checked.
      !! 
      !! The output is
      !! \[
      !!     \mathbf{U}^T \mathbf{V}, \qquad \text{and }  \qquad \mathbf{V_perp}^T \mathbf{V}
      !!     \hat{\mathbf{D}} = \mathbf{D} .
      !! \]
      real(dp),                   allocatable, intent(out) :: UTV(:,:)
      real(dp),                   allocatable, intent(out) :: VpTV(:,:)
      class(abstract_vector_rdp),              intent(in)  :: U(:)
      class(abstract_vector_rdp),              intent(in)  :: V(:)

      ! internals
      class(abstract_vector_rdp),             allocatable  :: Vp(:)
      real(dp),                               allocatable  :: wrk(:,:)
      integer :: ru, rv, r, info

      ru = size(U)
      rv = size(V)
      r  = ru + rv

      allocate(Vp(rv), source=V) ! Vp = V
      allocate(UTV( ru,rv)); UTV  = 0.0_dp
      allocate(VpTV(rv,rv)); VpTV = 0.0_dp

      ! orthonormalize second basis against first
      call orthogonalize_against_basis(Vp, U, info, if_chk_orthonormal=.false., beta=UTV)
      call check_info(info, 'orthogonalize_against_basis', this_module, 'project_onto_common_basis_rdp')
      allocate(wrk(rv,rv)); wrk = 0.0_dp
      call qr(Vp, wrk, info)
      call check_info(info, 'qr', this_module, 'project_onto_common_basis_rdp')

      ! compute inner product between second basis and its orthonormalized version
      VpTV = innerprod(Vp, V)

   end subroutine project_onto_common_basis_rdp

   real(dp) function increment_norm(X, U_lag, S_lag, ifnorm) result(inc_norm)
      !! This function computes the norm of the solution increment in a cheap way avoiding the
      !! construction of the full low-rank solutions.
      class(abstract_sym_low_rank_state_rdp)  :: X
      !! Low rank solution of current solution
      class(abstract_vector_rdp)              :: U_lag(:)
      !! Low-rank basis of lagged solution
      real(dp)                                :: S_lag(:,:)
      !! Coefficients of lagged solution
      logical, optional, intent(in) :: ifnorm
      logical                       :: ifnorm_
      !! Normalize solution by vector size?

      ! internals
      real(dp), dimension(:,:),                 allocatable :: D, V1, V2
      real(dp), dimension(:),                   allocatable :: svals       
      integer :: rk, rl

      ifnorm_ = optval(ifnorm, .true.)

      rk  = X%rk
      rl = size(U_lag)

      ! compute common basis
      call project_onto_common_basis_rdp(V1, V2, U_lag, X%U(:rk))

      ! project second low-rank state onto common basis and construct difference
      allocate(D(rk+rl,rk+rl)); D = 0.0_dp
      D(    :rl   ,      :rl   ) = S_lag - matmul(V1, matmul(X%S(:rk,:rk), transpose(V1)))
      D(rl+1:rl+rk,      :rl   ) =       - matmul(V2, matmul(X%S(:rk,:rk), transpose(V1)))
      D(    :rl   ,  rl+1:rl+rk) =       - matmul(V1, matmul(X%S(:rk,:rk), transpose(V2)))
      D(rl+1:rl+rk,  rl+1:rl+rk) =       - matmul(V2, matmul(X%S(:rk,:rk), transpose(V2)))

      ! compute Frobenius norm of difference
      inc_norm = sqrt(sum(svdvals(D))**2)
      if (ifnorm_) inc_norm = inc_norm/X%U(1)%get_size()

      return
   end function increment_norm

   real(dp) function low_rank_CALE_residual_norm(X, A, B, ifnorm) result(residual_norm)
      class(abstract_sym_low_rank_state_rdp) :: X
      !! Low-Rank factors of the solution.
      class(abstract_linop_rdp)              :: A
      !! Linear operator
      class(abstract_vector_rdp)             :: B(:)
      !! Low-Rank inhomogeneity.
      logical, optional                      :: ifnorm
      logical                                :: ifnorm_
      !! Normalize the norm by the vector size?
      ! internals
      integer :: i, rk, rkb, n, info
      class(abstract_vector_rdp), allocatable :: Q(:)
      real(dp), dimension(:,:), allocatable :: R, R_shuffle, sqrt_S
      ! optional arguments
      ifnorm_ = optval(ifnorm, .true.)

      rk  = X%rk
      rkb = size(B)
      n   = 2*rk + rkb
      allocate(Q(n), source=B(1)); call zero_basis(Q)
      allocate(R(n,n), R_shuffle(n,n)); R = 0.0_dp; R_shuffle = 0.0_dp

      ! fill the basis
      allocate(sqrt_S(rk,rk)); sqrt_S = 0.0_dp
      call sqrtm(X%S(:rk,:rk), sqrt_S, info)
      call check_info(info, 'sqrtm', this_module, 'low_rank_CALE_residual_norm')
      block
         class(abstract_vector_rdp), allocatable :: Xwrk(:)
         call linear_combination(Xwrk, X%U(:rk), sqrt_S)
         call copy(Q(rk+1:2*rk), Xwrk)
      end block
      do i = 1, rk
         call A%matvec(Q(rk+i), Q(i))
      end do
      call copy(Q(2*rk+1:), B(:))

      call qr(Q, R, info)
      call check_info(info, 'qr', this_module, 'low_rank_CALE_residual_norm')

      R_shuffle(:,      :  rk) = R(:,  rk+1:2*rk)
      R_shuffle(:,  rk+1:2*rk) = R(:,      :  rk)
      R_shuffle(:,2*rk+1:    ) = R(:,2*rk+1:    )

      residual_norm = norm2(matmul(R_shuffle, transpose(R)))
      if (ifnorm_) residual_norm = residual_norm/B(1)%get_size()
      
      return
   end function low_rank_CALE_residual_norm

   real(dp) function coefficient_matrix_norm(X, ifnorm) result(norm)
      !! This function computes the Frobenius norm of a low-rank approximation via an SVD of the (small) coefficient matrix
      class(abstract_sym_low_rank_state_rdp), intent(in) :: X
      !! Low rank solution of which to compute the norm
      logical, optional, intent(in) :: ifnorm
      logical                       :: ifnorm_
      !! Normalize the norm by the vector size?
      ifnorm_ = optval(ifnorm, .true.)
      norm = sqrt(sum(svdvals(X%S(:X%rk,:X%rk))**2))
      if (ifnorm_) norm = norm/X%U(1)%get_size()
      return
   end function coefficient_matrix_norm

   subroutine find_rank(found, irk, svals, tol)
      logical, intent(out) :: found
      !! Result of the check.
      integer, intent(out) :: irk
      !! Index of the first singular value below tolerance
      real(dp), intent(in) :: svals(:)
      !! Singular values to search
      real(dp), intent(in) :: tol
      !! Tolerance for the smallest resolved singular value

      found = .false.
      tol_chk: do irk = 1, size(svals)
         if ( svals(irk) < tol ) then
            found = .true.
            exit tol_chk
         end if
      end do tol_chk
   end subroutine find_rank

   subroutine increase_rank(X)
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.

      ! internal
      character(len=*), parameter :: this_procedure = 'increase_rank'
      integer :: rk, rkmax, info
      character(len=128) :: msg

      rkmax = size(X%U)
      if (rk == rkmax) then ! cannot increase rank without reallocating X%U and X%S
         write(msg,'(A,I0,A,A)') 'Cannot increase rank, rkmax = ', rkmax, ' is reached. ', &
                  & 'Increase rkmax and restart!'
         call stop_error(msg, this_module, this_procedure)
      else
         
         X%rk = X%rk + 1
         rk = X%rk ! this is only to make the code more readable
         ! set coefficients to zero (for redundancy)
         X%S(:rk, rk) = 0.0_dp 
         X%S( rk,:rk) = 0.0_dp
         ! add random vector ...
         call X%U(rk)%rand(.false.)
         ! ... and orthonormalize
         call orthogonalize_against_basis(X%U(rk), X%U(:rk-1), info, if_chk_orthonormal=.false.)
         call check_info(info, 'orthogonalize_against_basis', this_module, this_procedure)
         call X%U(rk)%scal(1.0_dp / X%U(rk)%norm())

      end if

   end subroutine increase_rank

   subroutine decrease_rank(X, U, svals, rk)
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      real(dp), intent(in) :: U(:,:)
      !! Left singular vectors of X
      real(dp), intent(in) :: svals(:)
      !! Singular values of S
      integer, intent(in) :: rk
      !! Desired output rank

      ! internal
      character(len=*), parameter :: this_procedure = 'decrease_rank'
      character(len=128) :: msg

      ! sanity check
      if (rk > size(X%U)) then
         write(msg,'(A,I0,1X,I0)') 'Invalid rank input: ', rk, size(X%U)
         call stop_error(msg, this_module, this_procedure)
      end if

      ! rotate basis onto principal axes
      block
         class(abstract_vector_rdp), allocatable :: Xwrk(:)
         call linear_combination(Xwrk, X%U(:rk), U(:rk,:rk))
         call copy(X%U(:rk), Xwrk)
      end block
      X%S(:rk,:rk) = diag(svals(:rk))

   end subroutine decrease_rank

   subroutine print_svals(X, svals, svals_lag, istep, nsteps)
      class(abstract_sym_low_rank_state_rdp), intent(in) :: X
      !! Low-Rank factors of the solution.
      real(dp), dimension(:), intent(in) :: svals
      !! Current singular values
      real(dp), dimension(:), intent(in) :: svals_lag
      !! Lagged singular values
      integer, intent(in) :: istep
      !! Current step
      integer, intent(in) :: nsteps
      !! Total number of steps

      ! Internal variables
      integer, parameter                  :: iline = 4       ! # data points per line
      integer                             :: i, j, is, ie, irk, ifmt, irkfmt
      character(len=128)                  :: msg, fmt_sval
      real(dp), dimension(:), allocatable :: dsvals

      ! Pretty output
      ifmt = max(5,ceiling(log10(real(nsteps))))
      irkfmt = max(3,ceiling(log10(real(size(X%U)))))
      write(fmt_sval,'(A,5(I0,A))') '("Step ",I', ifmt, ',"/",I', ifmt, ',": T= ",F10.4,1X,I', irkfmt, &
               & '" : ",A,"[",I', irkfmt, ',"-",I', irkfmt, ',"]",*(E12.5))'
      
      irk = min(size(svals), size(svals_lag))
      allocate(dsvals(irk)); dsvals = 0.0_dp
      do i = 1, irk
         dsvals(i) = abs(svals(i)-svals_lag(i))/svals(i)
      end do
      do i = 1, ceiling(float(X%rk)/iline)
         is = (i-1)*iline+1; ie = min(X%rk,i*iline)
         write(msg,fmt_sval) istep, nsteps, X%tot_time, X%rk, " SVD abs", is, ie, ( svals(j), j = is, ie )
         call log_information(msg, this_module, 'DLRA_main')
      end do
      do i = 1, ceiling(float(irk)/iline)
         is = (i-1)*iline+1; ie = min(irk,i*iline)
         write(msg,fmt_sval) istep, nsteps, X%tot_time, X%rk, "dSVD rel", is, ie, ( dsvals(j) , j = is, ie )
         call log_information(msg, this_module, 'DLRA_main')
      end do
   end subroutine print_svals

   logical function is_converged(X, svals, svals_lag, opts, if_lastep) result(converged)
      !! This function checks the convergence of the solution based on the (relative) increment in the singular values
      class(abstract_sym_low_rank_state_rdp) :: X
      real(dp)                   :: svals(:)
      real(dp)                   :: svals_lag(:)
      type(dlra_opts)            :: opts
      logical                    :: if_lastep
      ! internals
      real(dp), allocatable :: dsvals(:)
      integer :: i
      real(dp) :: norm, norm_lag, dnorm
      character(len=128) :: msg, prefix
      character(len=128), parameter :: fmt = '(A,I8,F15.8,A,2(E15.7,1X),A,E15.7)'

      norm     = sqrt(sum(svals**2))
      norm_lag = sqrt(sum(svals_lag**2))

      allocate(dsvals(size(svals)))
      do i = 1, size(svals)
         dsvals(i) = abs(svals(i) - svals_lag(i))
      end do

      dnorm = sqrt(sum(dsvals**2))

      if (opts%relative_inc) dnorm = dnorm/norm

      converged = .false.

      if (mod(X%step,opts%chkstep) == 0) then
         write(msg,fmt) 'Check state: ', X%step, X%time, ' svals lag ', norm, norm_lag, ' inc_norm ', dnorm
         call log_message(msg, this_module, 'DLRA')
      else if (if_lastep) then
         write(msg,fmt) 'Final state: ', X%step, X%time, ' svals lag ', norm, norm_lag, ' inc_norm ', dnorm
         call log_message(msg, this_module, 'DLRA')
      end if
      if (dnorm < opts%inc_tol) converged = .true.

      return
   end function is_converged

end module LightROM_Utils
