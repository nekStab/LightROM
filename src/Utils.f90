module LightROM_Utils
   ! stdlib
   use stdlib_strings, only: padl
   use stdlib_linalg, only : eye, diag, svd, svdvals, is_symmetric
   use stdlib_optval, only : optval
   ! LightKrylov for Linear Algebra
   use LightKrylov
   use LightKrylov, only : dp, wp => dp
   use LightKrylov_Logger
   use LightKrylov_AbstractVectors
   use LightKrylov_BaseKrylov, only : orthogonalize_against_basis
   use LightKrylov_Utils, only : abstract_opts, sqrtm
   ! LightROM
   use LightROM_AbstractLTIsystems
   
   implicit none 

   private :: this_module
   character(len=*), parameter :: this_module     = 'LR_Utils'
   character(len=*), parameter :: logfile_SVD_abs = 'Lyap_SVD_abs.dat'
   character(len=*), parameter :: logfile_SVD_rel = 'Lyap_SVD_rel.dat'

   public :: dlra_opts
   public :: coefficient_matrix_norm, increment_norm, low_rank_CALE_residual_norm
   public :: is_converged
   public :: write_logfile_headers, stamp_logfiles, log_settings
   public :: project_onto_common_basis
   public :: Balancing_Transformation
   public :: ROM_Petrov_Galerkin_Projection
   public :: ROM_Galerkin_Projection

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

   type, extends(abstract_opts), public :: dlra_opts
      !! Options container for the (rank-adaptive) projector-splitting dynalical low-rank approximation
      !! integrator
      integer :: mode = 1
      !! Time integration mode. Only 1st order (Lie splitting - mode 1) and 
      !! 2nd order (Strang splitting - mode 2) are implemented. (default: 1)
      !
      ! CONVERGENCE CHECK
      !
      integer :: chkstep = 10
      !! Time step interval at which convergence is checked and runtime information is printed (default: 10)
      real(wp) :: chktime = 1.0_wp
      !! Simulation time interval at which convergence is checked and runtime information is printed (default: 1.0)
      logical :: chkctrl_time = .true.
      !! Use time instead of timestep control (default: .true.)
      real(wp) :: inc_tol = 1e-6_wp
      !! Tolerance on the increment for convergence (default: 1e-6)
      logical :: relative_inc = .true.
      !! Tolerance control: Use relative values for convergence (default = .true.)
      !
      ! RANK-ADPATIVE SPECIFICS
      !
      logical :: if_rank_adaptive = .true.
      !! Allow rank-adaptivity
      real(wp) :: tol = 1e-6_wp
      !! Tolerance on the extra singular value to determine rank-adaptation
      logical :: use_err_est = .false.
      !! Choose whether to base the tolerance on 'tol' or on the splitting error estimate
      integer :: err_est_step = 10
      !! Time step interval for recomputing the splitting error estimate (only of use_err_est = .true.)
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
         call copy(T(1:rkmin), Xwrk)
         call linear_combination(Xwrk, Xc(1:rkmin), matmul(V(1:rkmin,1:rkmin), diag(Sigma)))
         call copy(Tinv(1:rkmin), Xwrk)
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
      class(abstract_lti_system_rdp),   intent(inout)  :: LTI
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
      class(abstract_lti_system_rdp),   intent(inout)  :: LTI
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
      real(wp),                               allocatable  :: wrk(:,:)
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
      allocate(wrk(rv,rv)); wrk = 0.0_wp
      call qr(Vp, wrk, info)
      call check_info(info, 'qr', module=this_module, procedure='project_onto_common_basis_rdp')

      ! compute inner product between second basis and its orthonormalized version
      call innerprod(VpTV, Vp, V)

      return
   end subroutine project_onto_common_basis_rdp

   real(dp) function increment_norm(X, U_lag, S_lag, ifnorm) result(inc_norm)
      !! This function computes the norm of the solution increment in a cheap way avoiding the
      !! construction of the full low-rank solutions.
      class(abstract_sym_low_rank_state_rdp)  :: X
      !! Low rank solution of current solution
      class(abstract_vector_rdp)              :: U_lag(:)
      !! Low-rank basis of lagged solution
      real(wp)                                :: S_lag(:,:)
      !! Coefficients of lagged solution
      logical, optional, intent(in) :: ifnorm
      logical                       :: ifnorm_
      !! Normalize solution by vector size?

      ! internals
      real(wp), dimension(:,:),                 allocatable :: D, V1, V2
      real(wp), dimension(:),                   allocatable :: svals       
      integer :: rk, rl

      ifnorm_ = optval(ifnorm, .true.)

      rk  = X%rk
      rl = size(U_lag)

      ! compute common basis
      call project_onto_common_basis_rdp(V1, V2, U_lag, X%U(:rk))

      ! project second low-rank state onto common basis and construct difference
      allocate(D(rk+rl,rk+rl)); D = 0.0_wp
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
      call check_info(info, 'sqrtm', module=this_module, procedure='low_rank_CALE_residual_norm')
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
      call check_info(info, 'qr', module=this_module, procedure='low_rank_CALE_residual_norm')

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

   logical function is_converged(X, svals, svals_lag, opts, if_lastep) result(converged)
      !! This function checks the convergence of the solution based on the (relative) increment in the singular values
      class(abstract_sym_low_rank_state_rdp) :: X
      real(wp)                   :: svals(:)
      real(wp)                   :: svals_lag(:)
      type(dlra_opts)            :: opts
      logical                    :: if_lastep
      ! internals
      real(wp),      allocatable :: dsvals(:)
      integer :: i
      real(wp) :: norm, norm_lag, dnorm
      character(len=128) :: msg, prefix

      norm     = sqrt(sum(svals**2))
      norm_lag = sqrt(sum(svals_lag**2))

      allocate(dsvals(size(svals)))
      do i = 1, size(svals)
         dsvals(i) = abs(svals(i) - svals_lag(i))
      end do

      dnorm    = sqrt(sum(dsvals**2))

      if (opts%relative_inc) dnorm = dnorm/norm

      converged = .false.

      if (if_lastep) then
         prefix = 'DLRA check (final time)'
      else
         prefix = 'DLRA check'
      end if

      write(msg,'(I8, F15.8, A,2(E15.7,1X),A,E15.7)') X%step, X%time, ' svals lag ', norm, norm_lag, ' inc_norm ', dnorm
      call logger%log_message(msg, module=this_module, procedure=prefix)
      if (dnorm < opts%inc_tol) converged = .true.

      return
   end function is_converged

   subroutine check_options(chkstep, tau, X, opts)
      integer,                                 intent(out)   :: chkstep
      real(wp),                                intent(in)    :: tau
      class(abstract_sym_low_rank_state_rdp),  intent(inout) :: X
      type(dlra_opts),                         intent(inout) :: opts

      ! internal
      character(len=128) :: msg
      type(dlra_opts) :: opts_default
      opts_default = dlra_opts()
      !
      ! CONVERGENCE CHECK
      !
      if (opts%chkctrl_time) then
         if (opts%chktime <= 0.0_wp) then
            opts%chktime = opts_default%chktime
            write(msg,'(A,E12.5,A)') 'Invalid chktime. Reset to default (',  opts%chktime,')'
            call logger%log_warning(msg, module=this_module, procedure='DLRA_check_options')
         end if
         chkstep = max(1, NINT(opts%chktime/tau))
         write(msg,'(A,E12.5,A,I0,A)') 'Convergence check every ', opts%chktime, ' time units (', chkstep, ' steps)'
         call logger%log_information(msg, module=this_module, procedure='DLRA_check_options')
      else
         if (opts%chkstep <= 0) then
            opts%chkstep = opts_default%chkstep
            write(msg,'(A,I0,A)') "Invalid chktime. Reset to default (",  opts%chkstep,")"
            call logger%log_warning(msg, module=this_module, procedure='DLRA_check_options')
         end if
         chkstep = opts%chkstep
         write(msg,'(A,I0,A)') 'Convergence check every ', opts%chkstep, ' steps (based on steps).'
         call logger%log_information(msg, module=this_module, procedure='DLRA_check_options')
      end if
      return
   end subroutine check_options

   subroutine write_logfile_headers(n0)
      integer, intent(in) :: n0
      ! internals
      integer :: i
      ! SVD absolute
      open (1234, file=logfile_SVD_abs, status='replace', action='write')
      write (1234, '(A8,2(A15,1X),A4)', ADVANCE='NO') 'istep', 'time', 'lag', 'rk'
      do i = 1, n0
         write (1234, '(A13,I2.2,1X)', ADVANCE='NO') 's', i
      end do
      write (1234, *) ''; close (1234)
      ! dSVD relative
      open (1234, file=logfile_SVD_rel, status='replace', action='write')
      write (1234, '(A8,2(A15,1X),A4)', ADVANCE='NO') 'istep', 'time', 'lag', 'rk'
      do i = 1, n0
         write (1234, '(A13,I2.2,1X)', ADVANCE='NO') 'ds', i
      end do
      write (1234, *) ''; close (1234)
      return
   end subroutine write_logfile_headers

   subroutine stamp_logfiles(X, lag, svals, dsvals)
      class(abstract_sym_low_rank_state_rdp),  intent(in) :: X
      real(dp), intent(in) :: lag
      real(dp), dimension(:), intent(in) :: svals
      real(dp), dimension(:), intent(in) :: dsvals
      ! SVD absolute
      open (1234, file=logfile_SVD_abs, status='old', action='write', position='append')
      write (1234, '(I8,2(1X,F15.9),I4)', ADVANCE='NO') X%step, X%time, lag, X%rk
      write (1234, '(*(1X,F15.9))') svals
      close (1234)
      ! dSVD relative
      open (1234, file=logfile_SVD_rel, status='old', action='write', position='append')
      write (1234, '(I8,2(1X,F15.9),I4)', ADVANCE='NO') X%step, X%time, lag, X%rk
      write (1234, '(*(1X,F15.9))') dsvals
      close (1234)
      return
   end subroutine stamp_logfiles

   subroutine log_settings(X, Tend, tau, nsteps, opts)
      class(abstract_sym_low_rank_state_rdp),  intent(in) :: X
      real(dp), intent(in) :: Tend
      real(dp), intent(in) :: tau
      integer, intent(in) :: nsteps
      type(dlra_opts), intent(in) :: opts
      ! internals
      character(len=128) :: msg, ctype
      call logger%log_message('###### solver settings ######', module=this_module, procedure='DLRA')
      write(msg,'(A15," : ", F15.8)') padl('t0',15), X%time
      call logger%log_message(msg, module=this_module, procedure='DLRA')
      write(msg,'(A15," : ", F15.8)') padl('tf',15), Tend
      call logger%log_message(msg, module=this_module, procedure='DLRA')
      write(msg,'(A15," : ", F15.8)') padl('dt',15), tau
      call logger%log_message(msg, module=this_module, procedure='DLRA')
      write(msg,'(A15," : ", I0)')    padl('nsteps',15),  nsteps
      call logger%log_message(msg, module=this_module, procedure='DLRA')
      write(msg,'(A15," : ", I0)')    padl('t-order',15), opts%mode
      call logger%log_message(msg, module=this_module, procedure='DLRA')
      write(msg,'(A15," : ", L)')     padl('adaptive rank',15), opts%if_rank_adaptive
      call logger%log_message(msg, module=this_module, procedure='DLRA')
      if (opts%if_rank_adaptive) then
         write(msg,'(A15," : ", I0)') padl('rk_init',15), X%rk
         call logger%log_message(msg, module=this_module, procedure='DLRA')
         write(msg,'(A15," : ", I0)') padl('rk_max',15), size(X%U)
         call logger%log_message(msg, module=this_module, procedure='DLRA')
         write(msg,'(A15," : sigma_{r+1} < ", E15.8)') padl('adapt. tol.',15), opts%tol
         call logger%log_message(msg, module=this_module, procedure='DLRA')
      else
         write(msg,'(A15," : ", I0)') padl('rk',15), X%rk
         call logger%log_message(msg, module=this_module, procedure='DLRA')
      end if
      if (opts%relative_inc) then
         ctype = 'relative'
      else
         ctype = 'absolute'
      end if
      write(msg,'(A15," : ",A,A)')    padl('convergence',15), trim(ctype), ' increment of the solution 2-norm'
      call logger%log_message(msg, module=this_module, procedure='DLRA')
      write(msg,'(A15," : ", E15.8)') padl('tol',15), opts%inc_tol
      call logger%log_message(msg, module=this_module, procedure='DLRA')
      if (opts%chkctrl_time) then
         write(msg,'("  Output every ",F8.4," time units (",I0," steps)")') opts%chktime, int(opts%chktime/tau)
         call logger%log_message(msg, module=this_module, procedure='DLRA')
      else
         write(msg,'("  Output every ",I0," steps (",F8.4," time units)")') opts%chkstep, opts%chkstep*tau
         call logger%log_message(msg, module=this_module, procedure='DLRA')
      end if
      call logger%log_message('###### solver settings ######', module=this_module, procedure='DLRA')
   end subroutine

end module LightROM_Utils
