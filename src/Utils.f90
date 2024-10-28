module LightROM_Utils
   ! stdlib
   use stdlib_linalg, only : eye, diag, svd, svdvals, is_symmetric
   use stdlib_optval, only : optval
   ! LightKrylov for Linear Algebra
   use LightKrylov
   use LightKrylov, only : dp, wp => dp
   use LightKrylov_Logger
   use LightKrylov_AbstractVectors
   use LightKrylov_BaseKrylov, only : orthogonalize_against_basis
   use LightKrylov_Utils, only : abstract_opts
   ! LightROM
   use LightROM_AbstractLTIsystems
   
   implicit none 

   private :: this_module
   character(len=*), parameter :: this_module = 'LightROM_Utils'

   public :: dlra_opts, abstract_outpost_rdp
   public :: compute_norm, compute_increment_norm
   public :: is_converged
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
      !! Tolerance on the increment norm for convergence (default: 1e-6)
      logical :: relative_norm = .true.
      !! Tolerance control: Check convergence for dX/X (true) or dX (false)? (default: .true.)
      !
      ! OUTPUT CALLBACK
      !
      logical :: ifIO = .false.
      !! Oupost during simulation (default = .false.). The final result is always returned.
      integer :: iostep = 100
      !! Time step interval at which to outpost
      real(wp) :: iotime = 1.0_wp
      !! Simulation time interval at which convergence is checked and runtime information is printed (default: 1.0)
      logical :: ioctrl_time = .true.
      !! IO control: use time instead of timestep control (default: .true.)
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

   abstract interface
      subroutine abstract_outpost_rdp(X, info)
         import abstract_low_rank_representation_rdp
         !! Abstract interface to define the matrix exponential-vector product.
         class(abstract_low_rank_representation_rdp), intent(in)  :: X
         integer,                                     intent(out) :: info
      end subroutine abstract_outpost_rdp
   end interface

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

   subroutine compute_increment_norm(nrm, U, S, U_lag, S_lag, ifnorm)
      !! This function computes the norm of the solution increment in a cheap way avoiding the
      !! construction of the full low-rank solutions.
      real(wp),                               intent(out)   :: nrm
      !! Increment norm of current timestep
      class(abstract_vector_rdp),             intent(in)    :: U(:)
      !! Low-rank basis of current solution
      real(wp),                               intent(in)    :: S(:,:)
      !! Coefficients of current solution
      class(abstract_vector_rdp),             intent(in)    :: U_lag(:)
      !! Low-rank basis of lagged solution
      real(wp),                               intent(in)    :: S_lag(:,:)
      !! Coefficients of lagged solution
      logical, optional, intent(in) :: ifnorm
      logical                       :: ifnorm_
      !! Normalize solution by vector size?

      ! internals
      real(wp), dimension(:,:),                 allocatable :: D, V1, V2
      real(wp), dimension(:),                   allocatable :: svals       
      integer :: r, rl

      ifnorm_ = optval(ifnorm, .true.)

      r  = size(U)
      rl = size(U_lag)

      ! compute common basis
      call project_onto_common_basis_rdp(V1, V2, U_lag, U)

      ! project second low-rank state onto common basis and construct difference
      allocate(D(r+rl,r+rl)); D = 0.0_wp
      D(    :rl  ,     :rl  ) = S_lag - matmul(V1, matmul(S, transpose(V1)))
      D(rl+1:rl+r,     :rl  ) =       - matmul(V2, matmul(S, transpose(V1)))
      D(    :rl  , rl+1:rl+r) =       - matmul(V1, matmul(S, transpose(V2)))
      D(rl+1:rl+r, rl+1:rl+r) =       - matmul(V2, matmul(S, transpose(V2)))

      ! compute Frobenius norm of difference
      nrm = sqrt(sum(svdvals(D))**2)
      if (ifnorm_) nrm = nrm/U(1)%get_size()

      return
   end subroutine compute_increment_norm

   real(dp) function compute_norm(X, ifnorm) result(nrm)
      !! This function computes the Frobenius norm of a low-rank approximation via an SVD of the (small) coefficient matrix
      class(abstract_sym_low_rank_state_rdp), intent(in) :: X
      !! Low rank solution of which to compute the norm
      logical, optional, intent(in) :: ifnorm
      logical                       :: ifnorm_
      !! Normalize the norm by the vector size?
      ifnorm_ = optval(ifnorm, .true.)
      nrm = sqrt(sum(svdvals(X%S(:X%rk,:X%rk))**2))
      if (ifnorm_) nrm = nrm/X%U(1)%get_size()
   end function compute_norm

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

   subroutine get_step(chkstep, iostep, opts, tau)
      integer,         intent(out)   :: chkstep 
      integer,         intent(out)   :: iostep 
      type(dlra_opts), intent(inout) :: opts
      real(wp),        intent(in)    :: tau

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
            call logger%log_warning(msg, module=this_module, procedure='DLRA_get_step')
         end if
         chkstep = max(1, NINT(opts%chktime/tau))
         write(msg,'(A,E12.5,A,I0,A)') 'Convergence check every ', opts%chktime, ' time units (', chkstep, ' steps)'
         call logger%log_information(msg, module=this_module, procedure='DLRA_get_step')
      else
         if (opts%chkstep <= 0) then
            opts%chkstep = opts_default%chkstep
            write(msg,'(A,I0,A)') "Invalid chktime. Reset to default (",  opts%chkstep,")"
            call logger%log_warning(msg, module=this_module, procedure='DLRA_get_step')
         end if
         chkstep = opts%chkstep
         write(msg,'(A,I0,A)') 'Convergence check every ', chkstep, ' steps (based on steps).'
         call logger%log_information(msg, module=this_module, procedure='DLRA_get_step')
      end if
      !
      ! OUTPUT CALLBACK
      !
      if (opts%ifIO) then ! callback activated
         if (opts%ioctrl_time) then
            if (opts%iotime <= 0.0_wp) then
               opts%iotime = opts_default%iotime
               write(msg,'(A,E12.5,A)') 'Invalid iotime. Reset to default (',  opts%iotime,')'
               call logger%log_warning(msg, module=this_module, procedure='DLRA_get_step')
            end if
            iostep = max(1, NINT(opts%iotime/tau))
            write(msg,'(A,E12.5,A,I0,A)') 'Output every ', opts%iotime, ' time units (', iostep, ' steps)'
            call logger%log_information(msg, module=this_module, procedure='DLRA_get_step')
         else
            if (opts%iostep <= 0) then
               opts%iostep = opts_default%iostep
               write(msg,'(A,I0,A)') "Invalid iotime. Reset to default (",  opts%iostep,")"
               call logger%log_warning(msg, module=this_module, procedure='DLRA_get_step')
            end if
            iostep = opts%iostep
            write(msg,'(A,I0,A)') 'Output every ', iostep, ' steps (based on steps).'
            call logger%log_information(msg, module=this_module, procedure='DLRA_get_step')
         end if
      else
         iostep = 0
         call logger%log_information('No runtime output.', module=this_module, procedure='DLRA_get_step')
      end if
      return
   end subroutine get_step

end module LightROM_Utils
