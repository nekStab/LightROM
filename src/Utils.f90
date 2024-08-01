module LightROM_Utils
   ! stdlib
   use stdlib_linalg, only : eye, diag, svd, svdvals, is_symmetric
   use stdlib_optval, only : optval
   use stdlib_stats_distribution_normal, only: normal => rvs_normal
   use stdlib_logger, only : logger => global_logger
   ! LightKrylov for Linear Algebra
   use LightKrylov
   use LightKrylov_Constants
   use LightKrylov, only : dp, wp => dp
   use LightKrylov_Logger
   use LightKrylov_AbstractVectors
   use LightKrylov_BaseKrylov, only : orthonormalize_basis, orthogonalize_against_basis
   use LightKrylov_Utils, only : abstract_opts, assert_shape, sqrtm, eigh
   ! LightROM
   use LightROM_AbstractLTIsystems
   
   implicit none 

   private :: this_module
   character(len=*), parameter :: this_module = 'LightROM_Utils'

   ! Balancing Transformations and projections   
   public :: Balancing_Transformation
   public :: ROM_Petrov_Galerkin_Projection
   public :: ROM_Galerkin_Projection
   public :: project_onto_common_basis

   ! Utilities for matrix norm computations
   public :: dense_frobenius_norm
   public :: increment_norm
   public :: CALE_res_norm

   ! Miscellenous utils
   public :: chk_opts
   public :: is_converged
   public :: random_low_rank_state

   interface Balancing_Transformation
      module procedure Balancing_Transformation_rdp
   end interface

   interface ROM_Petrov_Galerkin_Projection
      module procedure ROM_Petrov_Galerkin_Projection_rdp
   end interface

   interface ROM_Galerkin_Projection
      module procedure ROM_Galerkin_Projection_rdp
   end interface

   interface project_onto_common_basis
      module procedure project_onto_common_basis_rdp
   end interface

   interface random_low_rank_state
      module procedure random_low_rank_state_rdp
   end interface

   type, extends(abstract_opts), public :: dlra_opts
      !! Options container for the (rank-adaptive) projector-splitting dynalical low-rank approximation
      !! integrator
      integer :: mode = 1
      !! Time integration mode. Only 1st order (Lie splitting - mode 1) and 
      !! 2nd order (Strang splitting - mode 2) are implemented. (default: 1)
      logical :: verbose = .false.
      !! Verbosity control (default: .false.)
      !
      !! CHKSTEP
      integer :: chkstep = 10
      !! Time step interval at which convergence is checked and runtime information is printed (default: 10)
      real(wp) :: chktime = 1.0_wp
      !! Simulation time interval at which convergence is checked and runtime information is printed (default: 1.0)
      logical :: chkctrl_time = .true.
      !! IO control: use time instead of timestep control (default: .true.)
      logical :: print_svals = .true.
      !! IO control: compute and print singular values of solution (default: .true.)
      !
      !! INCREMENT NORM
      logical :: chk_convergence = .true.
      !! Toggle whether to check for convergence of the solution
      real(wp) :: inc_tol = 1e-6_wp
      !! Tolerance on the increment norm for convergence (default: 1e-6)
      logical :: relative_norm = .true.
      !! Tolerance control: Check convergence for dX/X (true) or dX (false)? (default: .true.)
      logical :: if_rank_adaptive = .true.
      !! Allow rank-adaptivity
      !
      ! RANK-ADPATIVE SPECIFICS
      !    
      !! INITIAL RANK
      integer :: ninit = 10
      !! Number of time steps to run the integrator when determining the initial rank (default: 10)
      real(wp) :: tinit = 1.0_wp
      !! Physical time to run the integrator when determining the initial rank (default: 1.0)
      logical :: initctrl_step = .true.
      !! Init control: use ninit to determine the integration time for initial rank (default: .true.)
      !
      !! TOLERANCE
      real(wp) :: tol = 1e-6_wp
      !! Tolerance on the extra singular value to determine rank-adaptation
      integer :: err_est_step = 10
      !! Time step interval for recomputing the splitting error estimate (only of use_err_est = .true.)
      logical :: use_err_est = .false.
      !! Choose whether to base the tolerance on 'tol' or on the splitting error estimate
   end type

contains

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
      real(wp),                            intent(out)   :: S(:)
      !! Singular values of the BT
      class(abstract_vector_rdp),          intent(out)   :: Tinv(:)
      !! Inverse balancing transformation
      class(abstract_vector_rdp),          intent(in)    :: Xc(:)
      !! Low-rank representation of the Controllability Gramian
      class(abstract_vector_rdp),          intent(in)    :: Yo(:)
      !! Low-rank representation of the Observability Gramian

      ! internal variables
      integer                                :: i, rkc, rko, rk, rkmin
      real(wp),                  allocatable :: LRCrossGramian(:,:)
      real(wp),                  allocatable :: Swrk(:,:)
      real(wp),                  allocatable :: Sigma(:)
      real(wp),                  allocatable :: V(:,:), W(:,:)

      rkc   = size(Xc)
      rko   = size(Yo)
      rk    = max(rkc, rko)
      rkmin = min(rkc, rko) 

      ! compute inner product with Gramian bases and compte SVD
      allocate(LRCrossGramian(rkc,rko)); allocate(V(rko,rko)); allocate(W(rkc,rkc))
      call innerprod(LRCrossGramian, Xc, Yo)
      call svd(LRCrossGramian, S, V, W)

      allocate(Sigma(rkmin))
      do i = 1, rkmin
         Sigma(i) = 1/sqrt(S(i))
      enddo
      block
         class(abstract_vector_rdp), allocatable :: Xwrk(:)
         call linear_combination(Xwrk, Yo(1:rkmin), matmul(W(1:rkmin,1:rkmin), diag(Sigma)))
         call copy_basis(T(1:rkmin), Xwrk)
         call linear_combination(Xwrk, Xc(1:rkmin), matmul(V(1:rkmin,1:rkmin), diag(Sigma)))
         call copy_basis(Tinv(1:rkmin), Xwrk)
      end block
         
      return
   end subroutine Balancing_Transformation_rdp

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
      real(wp),            allocatable, intent(out)    :: Ahat(:, :)
      !! Reduced-order dynamics matrix.
      real(wp),            allocatable, intent(out)    :: Bhat(:, :)
      !! Reduced-order input-to-state matrix.
      real(wp),            allocatable, intent(out)    :: Chat(:, :)
      !! Reduced-order state-to-output matrix.
      real(wp),            allocatable, intent(out)    :: D(:, :)
      !! Feed-through matrix
      class(abstract_lti_system_rdp),   intent(in)     :: LTI
      !! Large-scale LTI to project
      class(abstract_vector_rdp),       intent(in)     :: T(:)
      !! Balancing transformation
      class(abstract_vector_rdp),       intent(in)     :: Tinv(:)
      !! Inverse balancing transformation

      ! internal variables
      integer                                          :: i, rk, rkc, rkb
      class(abstract_vector_rdp),       allocatable    :: Uwrk(:)
      real(wp),                         allocatable    :: Cwrk(:, :)

      rk  = size(T)
      rkb = size(LTI%B)
      rkc = size(LTI%CT)
      allocate(Uwrk(rk), source=T(1)); call zero_basis(Uwrk)
      allocate(Ahat(1:rk, 1:rk ));                  Ahat = 0.0_wp
      allocate(Bhat(1:rk, 1:rkb));                  Bhat = 0.0_wp
      allocate(Cwrk(1:rk, 1:rkc));                  Cwrk = 0.0_wp
      allocate(Chat(1:rkc,1:rk ));                  Chat = 0.0_wp
      allocate(D(1:size(LTI%D,1),1:size(LTI%D,2))); D    = 0.0_wp

      do i = 1, rk
         call LTI%A%matvec(Tinv(i), Uwrk(i))
      end do
      call innerprod(Ahat, T, Uwrk)
      call innerprod(Bhat, T, LTI%B)
      call innerprod(Cwrk, LTI%CT, Tinv)
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
      real(wp),            allocatable, intent(out)    :: Ahat(:, :)
      !! Reduced-order dynamics matrix.
      real(wp),            allocatable, intent(out)    :: Bhat(:, :)
      !! Reduced-order input-to-state matrix.
      real(wp),            allocatable, intent(out)    :: Chat(:, :)
      !! Reduced-order state-to-output matrix.
      real(wp),            allocatable, intent(out)    :: D(:, :)
      !! Feed-through matrix
      class(abstract_lti_system_rdp),   intent(in)     :: LTI
      !! Large-scale LTI to project
      class(abstract_vector_rdp),       intent(inout)  :: T(:)
      !! Balancing transformation

      call ROM_Petrov_Galerkin_Projection(Ahat, Bhat, Chat, D, LTI, T, T)

      return
   end subroutine ROM_Galerkin_Projection_rdp

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
      real(wp),                   allocatable, intent(out) :: UTV(:,:)
      real(wp),                   allocatable, intent(out) :: VpTV(:,:)
      class(abstract_vector_rdp),              intent(in)  :: U(:)
      class(abstract_vector_rdp),              intent(in)  :: V(:)

      ! internals
      class(abstract_vector_rdp),             allocatable  :: Vp(:)
      integer :: ru, rv, r, info

      ru = size(U)
      rv = size(V)
      r  = ru + rv

      allocate(Vp(rv), source=V) ! Vp = V
      allocate(UTV( ru,rv)); UTV  = 0.0_wp
      allocate(VpTV(rv,rv)); VpTV = 0.0_wp

      ! orthonormalize second basis against first
      call orthogonalize_against_basis(Vp, U, info, if_chk_orthonormal=.false., beta=UTV)
      call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='project_onto_common_basis_rdp')
      call orthonormalize_basis(Vp)
      call check_info(info, 'qr', module=this_module, procedure='project_onto_common_basis_rdp')

      ! compute inner product between second basis and its orthonormalized version
      call innerprod(VpTV, Vp, V)

      return
   end subroutine project_onto_common_basis_rdp

   real(wp) function dense_frobenius_norm(S, scale) result(nrm)
      !! This function computes the Frobenius norm of a low-rank approximation via an SVD of the (small) coefficient matrix
      real(wp), intent(in) :: S(:,:)
      !! Dense (coefficient) matrix to compute frobenius norm of
      real(wp), optional, intent(in) :: scale
      real(wp)                       :: scale_
      !! Scaling factor
      
      scale_ = optval(scale, 1.0_wp)
      if (scale_ < atol_dp) call stop_error('Wrong input for scale', module=this_module, procedure='compute_norm')
      
      nrm = sqrt(sum(svdvals(S)**2))/scale_

   end function dense_frobenius_norm

   real(wp) function increment_norm(U, S, U_lag, S_lag, scale) result(nrm)
      !! This function computes the norm of the solution increment in a cheap way avoiding the
      !! construction of the full low-rank solutions.
      class(abstract_vector_rdp),             intent(in)    :: U(:)
      !! Low-rank basis of current solution
      real(wp),                               intent(in)    :: S(:,:)
      !! Coefficients of current solution
      class(abstract_vector_rdp),             intent(in)    :: U_lag(:)
      !! Low-rank basis of lagged solution
      real(wp),                               intent(in)    :: S_lag(:,:)
      !! Coefficients of lagged solution
      real(wp),                     optional, intent(in)    :: scale
      real(wp)                                              :: scale_
      !! Scaling factor

      ! internals
      real(wp), dimension(:,:),                 allocatable :: D, V1, V2
      integer :: r, rl

      scale_ = optval(scale, 1.0_wp)
      if (scale_ < atol_dp) call stop_error('Wrong input for scale', module=this_module, procedure='compute_increment_norm')

      r  = size(U)
      rl = size(U_lag)

      ! compute common basis
      call project_onto_common_basis(V1, V2, U_lag, U)

      ! project second low-rank state onto common basis and construct difference
      allocate(D(r+rl,r+rl)); D = 0.0_wp
      D(    :rl  ,     :rl  ) = S_lag - matmul(V1, matmul(S, transpose(V1)))
      D(rl+1:rl+r,     :rl  ) =       - matmul(V2, matmul(S, transpose(V1)))
      D(    :rl  , rl+1:rl+r) =       - matmul(V1, matmul(S, transpose(V2)))
      D(rl+1:rl+r, rl+1:rl+r) =       - matmul(V2, matmul(S, transpose(V2)))

      ! compute Frobenius norm of difference
      nrm = dense_frobenius_norm(D, scale_)

   end function increment_norm

   real(dp) function CALE_res_norm(X, A, B, scale) result(res)
      !! This function computes the Frobenius norm of a low-rank approximation via an SVD of the (small) coefficient matrix
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_linop_rdp),              intent(in) :: A
      !! Linear operator of the Lyapunov equation
      class(abstract_vector_rdp),             intent(in) :: B(:)
      !! Forcing term of the Lyapunov equation
      real(wp), optional, intent(in) :: scale
      real(wp)                       :: scale_
      !! Scaling factor
      
      ! internals
      integer :: i, info, rk
      class(abstract_vector_rdp), allocatable :: Q(:)   ! partial basis
      real(wp), dimension(:),    allocatable :: svals
      real(wp), dimension(:,:),  allocatable :: R, Rperm ! coefficient matrix of QR decomposition, sqrt(S)

      real(wp) :: V(X%rk,X%rk), lambda(X%rk)
      ! optional inputs
      scale_ = optval(scale, 1.0_wp)
      if (scale_ < atol_dp) call stop_error('Wrong input for scale', module=this_module, procedure='compute_norm')

      rk  = 2*X%rk + size(B)
      allocate(R(rk,rk), Rperm(rk,rk), svals(rk)); R = 0.0_wp; Rperm = 0.0_wp; svals = 0.0_wp
      allocate(Q(rk), source=X%U(1)); call zero_basis(Q)

      ! First X%rk columns are X%U @ sqrtm(X%S) 
      block
          class(abstract_vector_rdp), allocatable :: Xwrk(:)
          real(wp) :: Rwrk(X%rk,X%rk)
          call sqrtm(X%S(:X%rk,:X%rk), Rwrk, info)
          call check_info(info, 'sqrtm', module=this_module, procedure='compute_CALE_residual')
          call linear_combination(Xwrk, X%U(:X%rk), Rwrk)
          call copy_basis(Q(:X%rk), Xwrk)
      end block
      ! Second set of X%rk columns are given by A @ X%U @ sqrtm(X%S)
      do i = 1, X%rk
         call A%matvec(Q(i), Q(X%rk+i))
      end do
      ! Last columns are B
      call copy_basis(Q(2*X%rk+1:), B)
     
      ! compute QR decomposiion
      call qr(Q, R, info)
      call check_info(info, 'qr', module=this_module, procedure='compute_CALE_residual')
      
      ! Shuffle columns around
      Rperm(:,        :  X%rk) = R(:,  X%rk+1:2*X%rk)
      Rperm(:,  X%rk+1:2*X%rk) = R(:,        :  X%rk)
      Rperm(:,2*X%rk+1:      ) = R(:,2*X%rk+1:      ) ! last size(B) columns do not change
      
      ! compute norm of inner product
      R = matmul(Rperm, transpose(R))
      res = dense_frobenius_norm(R, scale_)

   end function CALE_res_norm

   logical function is_converged(nrm, nrmX, opts) result(converged)
      !! This function checks the convergence of the solution based on the (relative) increment norm
      real(wp),                   intent(in) :: nrm
      real(wp),         optional, intent(in) :: nrmX
      real(wp)                               :: nrmX_
      type(dlra_opts),  optional, intent(in) :: opts
      type(dlra_opts)                        :: opts_

      ! internals
      character*128 :: msg

      if (present(opts)) then
         opts_ = opts
      else
         opts_ = dlra_opts()
      end if

      if (present(nrmX)) then
         nrmX_ = nrmX
      else
         nrmX_ = 1.0_wp
      end if

      converged = .false.

      if (opts%relative_norm) then
         if (nrm/nrmX_ < opts%inc_tol) converged = .true.
      else
         if (nrm < opts%inc_tol) converged = .true.
      end if

   end function is_converged

   subroutine read_opts(opts, tau, rank_is_initialized)

      type(dlra_opts), intent(inout) :: opts
      real(wp),        intent(in)    :: tau
      logical,         intent(in)    :: rank_is_initialized

      ! internal
      character(len=128) :: msg
      type(dlra_opts) :: opts_default

      opts_default = dlra_opts()

      ! mode
      if ( opts%mode > 2 ) then
         opts%mode = 2
         write(msg, *) "Time-integration order for the operator splitting of d > 2 &
                      & requires adjoint solves and is not implemented. Resetting torder = 2." 
         if (io_rank()) call logger%log_warning(trim(msg), module=this_module, procedure='DLRA chk_opts')
      else if ( opts%mode < 1 ) then
         write(msg, '(A,I2)') "Invalid time-integration order specified: ", opts%mode
         call stop_error(trim(msg), module=this_module, procedure='DLRA chk_opts')
      endif

      ! chkctrl -- chkstep
      if (opts%chkctrl_time) then
         if (opts%chktime <= 0.0_wp) then
            opts%chktime = opts_default%chktime
            write(msg, '(A,F0.2,A,F0.2,A)') "Invalid chktime ( ", opts%chktime, " ). Reset to default ( ",  opts%chktime," )"
            if (io_rank()) call logger%log_warning(trim(msg), module=this_module, procedure='DLRA chk_opts')
         end if
         opts%chkstep = max(1, NINT(opts%chktime/tau))
         if (opts%verbose) then
            write(msg, '(A,F0.2,A,I6,A)') 'Output every ', opts%chktime, ' time units ( ', opts%chkstep, ' steps)'
            if (io_rank()) call logger%log_message(trim(msg), module=this_module, procedure='DLRA chk_opts')
         end if
      else
         if (opts%chkstep <= 0) then
            opts%chkstep = opts_default%chkstep
            write(msg, '(A,F0.2,A,I4,A)') "Invalid chktime ( ", opts%chktime, " ). Reset to default ( ",  opts%chkstep," )"
            if (io_rank()) call logger%log_message(trim(msg), module=this_module, procedure='DLRA chk_opts')
         end if
         if (opts%verbose) then
            write(msg,'(A,I6,A)') 'Output every ', opts%chkstep, ' steps (based on steps).'
            if (io_rank()) call logger%log_message(trim(msg), module=this_module, procedure='DLRA chk_opts')
         end if
      end if

      if (opts%if_rank_adaptive) then
         ! initctrl --> ninit
         if (.not.rank_is_initialized) then
            if (opts%verbose) then
               write(msg, '(A,E8.2)') 'Tolerance for rank adaptation: ', opts%tol
               if (io_rank()) call logger%log_message(trim(msg), module=this_module, procedure='DLRA chk_opts')
            end if
            if (opts%initctrl_step) then
               if (opts%ninit <= 0) then
                  opts%ninit = opts_default%ninit
                  write(msg, '(A,I4,A,I4,A)') "Invalid ninit ( ", opts%ninit, " ). Reset to default ( ",  opts%ninit," )"
                  if (io_rank()) call logger%log_warning(trim(msg), module=this_module, procedure='DLRA chk_opts')
               end if 
               if (opts%verbose) then
                  write(msg, '(A,I4,A)') 'Initial rank computed over ', opts%ninit, ' steps.'
                  if (io_rank()) call logger%log_message(trim(msg), module=this_module, procedure='DLRA chk_opts')
               end if
            else
               if (opts%tinit <= 0.0_wp) then
                  opts%tinit = opts_default%tinit
                  write(msg, '(A,F0.2,A,F0.2,A)') "Invalid tinit ( ", opts%tinit, " ). Reset to default ( ",  opts%tinit," )"
                  if (io_rank()) call logger%log_warning(trim(msg), module=this_module, procedure='DLRA chk_opts')
               end if
               opts%ninit = max(5, NINT(opts%tinit/tau))
               opts%tinit = opts%ninit*tau
               if (opts%verbose) then
                  write(msg, '(A,F0.2,A,I4,A)') 'Initial rank computed over ', opts%tinit, ' time units ( ', opts%ninit, ' steps)'
                  if (io_rank()) call logger%log_message(trim(msg), module=this_module, procedure='DLRA chk_opts')
               end if
            end if
         else
            if (opts%verbose) then
               write(msg, '(A)') 'Initial rank already set.'
               if (io_rank()) call logger%log_message(trim(msg), module=this_module, procedure='DLRA chk_opts')
            end if
         end if
      end if
      return
   end subroutine read_opts

   subroutine chk_opts(opts)

      type(dlra_opts), intent(inout) :: opts

      ! internal
      character(len=128) :: msg
      type(dlra_opts) :: opts_default

      opts_default = dlra_opts()

      ! mode
      if ( opts%mode > 2 ) then
         opts%mode = 2
         write(msg, *) "Time-integration order for the operator splitting of d > 2 &
                      & requires adjoint solves and is not implemented. Resetting torder = 2." 
         if (io_rank()) call logger%log_warning(trim(msg), module=this_module, procedure='DLRA chk_opts')
      else if ( opts%mode < 1 ) then
         write(msg, '(A,I2)') "Invalid time-integration order specified: ", opts%mode
         call stop_error(trim(msg), module=this_module, procedure='DLRA chk_opts')
      endif

      ! chkctrl -- chkstep
      if (opts%chkctrl_time) then
         if (opts%chktime <= 0.0_wp) then
            opts%chktime = opts_default%chktime
            write(msg, '(A,F0.2,A,F0.2,A)') "Invalid chktime ( ", opts%chktime, " ). Reset to default ( ",  opts%chktime," )"
            if (io_rank()) call logger%log_warning(trim(msg), module=this_module, procedure='DLRA chk_opts')
         end if
      else
         if (opts%chkstep <= 0) then
            opts%chkstep = opts_default%chkstep
            write(msg, '(A,F0.2,A,I4,A)') "Invalid chktime ( ", opts%chktime, " ). Reset to default ( ",  opts%chkstep," )"
            if (io_rank()) call logger%log_message(trim(msg), module=this_module, procedure='DLRA chk_opts')
         end if
      end if

      ! initctrl --> ninit
      if (opts%initctrl_step) then
         if (opts%ninit <= 0) then
            opts%ninit = opts_default%ninit
            write(msg, '(A,I4,A,I4,A)') "Invalid ninit ( ", opts%ninit, " ). Reset to default ( ",  opts%ninit," )"
            if (io_rank()) call logger%log_warning(trim(msg), module=this_module, procedure='DLRA chk_opts')
         end if 
      else
         if (opts%tinit <= 0.0_wp) then
            opts%tinit = opts_default%tinit
            write(msg, '(A,F0.2,A,F0.2,A)') "Invalid tinit ( ", opts%tinit, " ). Reset to default ( ",  opts%tinit," )"
            if (io_rank()) call logger%log_warning(trim(msg), module=this_module, procedure='DLRA chk_opts')
         end if
      end if
      return
   end subroutine chk_opts

   subroutine random_low_rank_state_rdp(U, S, V)
      class(abstract_vector_rdp),           intent(inout) :: U(:)
      real(dp),                             intent(inout) :: S(:,:)
      class(abstract_vector_rdp), optional, intent(inout) :: V(:)

      ! internals
      integer :: i, rk
      real(dp), dimension(:,:), allocatable :: mu, var

      rk = size(S, 1)
      call assert_shape(S, [rk, rk], 'random_low_rank_state', 'S')
      if (size(U) /= rk) call stop_error('Input basis U and coefficient matrix S have incompatible sizes', &
                                          & module=this_module, procedure='random_low_rank_state_rdp')

      allocate(mu(rk,rk), var(rk,rk))
      mu  = zero_rdp
      var = one_rdp
      S = normal(mu, var)
      
      call zero_basis(U)
      do i = 1, size(U)
         call U(i)%rand(.false.)
      end do
      call orthonormalize_basis(U)

      if (present(V)) then
         if (size(V) /= rk) call stop_error('Input basis V and coefficient matrix S have incompatible sizes', &
                                          & module=this_module, procedure='random_low_rank_state_rdp')
         call zero_basis(V)
         do i = 1, size(V)
            call V(i)%rand(.false.)
         end do                    
         call orthonormalize_basis(V)
      else
         ! symmetric
         S = 0.5*(S + transpose(S))
      end if

   end subroutine random_low_rank_state_rdp
   
end module LightROM_Utils
