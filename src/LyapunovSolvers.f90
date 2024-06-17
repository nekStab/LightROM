module LightROM_LyapunovSolvers
   !! This module provides the implementation of the Krylov-based solvers for the Differential Lyapunov
   !! equation based on the dynamic low-rank approximation and operator splitting.
   ! Standard library
   use stdlib_linalg, only : eye
   use stdlib_optval, only : optval
   ! LightKrylov modules
   use LightKrylov
   use LightKrylov, only: wp => dp
   use LightKrylov_Logger
   use LightKrylov_AbstractVectors
   use LightKrylov_ExpmLib
   use LightKrylov_BaseKrylov
   ! LightROM modules
   use LightROM_AbstractLTIsystems
   use LightROM_LyapunovUtils
   use LightROM_Utils
   
   implicit none

   ! global scratch arrays
   class(abstract_vector_rdp),  allocatable   :: U1(:)
   class(abstract_vector_rdp),  allocatable   :: Uwrk(:)
   class(abstract_vector_rdp),  allocatable   :: BBTU(:)
   real(wp),                    allocatable   :: Swrk(:,:)
   ! svd
   real(wp),                    allocatable   :: ssvd(:)
   real(wp),                    allocatable   :: Usvd(:,:), VTsvd(:,:)

   ! module name
   private :: this_module
   character*128, parameter :: this_module = 'LightKrylov_LyapunovSolvers'

   public :: numerical_low_rank_splitting_lyapunov_integrator
   public :: M_forward_map
   public :: G_forward_map_lyapunov
   public :: K_step_lyapunov
   public :: S_step_lyapunov
   public :: L_step_lyapunov

   interface numerical_low_rank_splitting_lyapunov_integrator
      module procedure numerical_low_rank_splitting_lyapunov_integrator_rdp
   end interface

   interface M_forward_map
      module procedure M_forward_map_rdp
   end interface

   interface G_forward_map_lyapunov
      module procedure G_forward_map_lyapunov_rdp
   end interface

   interface K_step_lyapunov
      module procedure K_step_lyapunov_rdp
   end interface

   interface S_step_lyapunov
      module procedure S_step_lyapunov_rdp
   end interface

   interface L_step_lyapunov
      module procedure L_step_lyapunov_rdp
   end interface

   contains

   subroutine numerical_low_rank_splitting_lyapunov_integrator_rdp(X, A, B, Tend, tau, torder, info, &
                                                                    & exptA, iftrans, ifverb, ifrk)
      !! Numerical integrator for the matrix-valued differential Lyapunov equation of the form
      !!
      !!    $$ \dot{\mathbf{X}} = \mathbf{A} \mathbf{X} + \mathbf{X} \mathbf{A}^T + \mathbf{B} \mathbf{B}^T $$
      !!
      !! where \( \mathbf{A} \) is a (n x n) Hurwitz matrix, \( \mathbf{X} \) is SPD and 
      !! \( \mathbf{B} \mathbf{B}^T \) is a rank-m rhs (m<<n).
      !!
      !! Since \( \mathbf{A} \) is Hurwitz, the equations converges to steady state for \( t \to \infty \), 
      !! which corresponds to the associated algebraic Lyapunov equation of the form
      !!
      !!    $$ \mathbf{0} = \mathbf{A} \mathbf{X} + \mathbf{X} \mathbf{A}^T + \mathbf{B} \mathbf{B}^T $$
      !!
      !! The algorithm is based on three main ideas:
      !!
      !! - The operator splitting scheme proposed by Lubich & Oseledets (2014) that splits the 
      !!   right-hand side of the differential equation into a linear stiff part that is solved
      !!   explicitly and a possibly non-linear non-stiff part which is solved numerically. The
      !!   two operators are then composed to obtain the integrator for the full Lyapunov equation.
      !! - The Dynamic Low-Rank Approximation for the solution of general matrix differential 
      !!   equations proposed by Nonnenmacher & Lubich (2007) which seeks to integrate only the
      !!   leading low-rank factors of the solution to a large system by updating the matrix 
      !!   factorization. The dynamical low-rank approximation scheme for the low-rank factors 
      !!   of the solution is itself solved using a projector-splitting technique to cheaply 
      !!   maintain orthonormality or the low-rank basis without explicit SVDs. 
      !! - This algorithm has been applied to the Lyapunov and Riccati equations by Mena et al. 
      !!   (2018) with improvements taking advantage of the symmetry of the problem/solution.
      !!
      !! **Algorithmic Features**
      !! 
      !! - Separate integration of the stiff inhomogeneous part of the Lyapunov equation and the
      !!   non-stiff inhomogeneity
      !! - Rank preserving time-integration that maintains orthonormality of the factorization
      !!   basis
      !! - The stiff part of the problem is solved using a time-stepper approach to approximate 
      !!   the action of the exponential propagator
      !!
      !! **Advantages**
      !!
      !! - Rank of the approximate solution is user defined
      !! - The timesteps of the stiff and non-stiff parts of the code are independent
      !! - The integrator is adjoint-free
      !! - The operator of the homogeneous part and the inhomogeneity are not needed explicitly
      !!   i.e. the algorithm is amenable to solution using Krylov methods (in particular for 
      !!   the solution of the stiff part of the problem)
      !! - No SVDs are necessary for this alogorithm
      !! - Lie and Strang splitting implemented allowing for first and second order integration
      !!   in time
      !!
      !! ** Limitations**
      !!
      !! - Rank of the approximate solution is user defined. The appropriateness of this 
      !!   approximation is not considered
      !! - The current implementation does not require an adjoint integrator. This means that
      !!   the temporal order of the basic operator splitting scheme is limited to 1 (Lie-Trotter
      !!   splitting) or at most 2 (Strang splitting). Higher order integrators are possible, but 
      !!   require at least some backward integration (via the adjoint) in BOTH parts of the splitting. 
      !!   (see Sheng-Suzuki and Goldman-Kaper theorems)
      !!
      !! **References**
      !! 
      !! - Koch, O.,Lubich, C. (2007). "Dynamical Low‐Rank Approximation", SIAM Journal on Matrix Analysis 
      !!   and Applications 29(2), 434-454
      !! - Lubich, C., Oseledets, I.V. (2014). "A projector-splitting integrator for dynamical low-rank 
      !!   approximation", BIT Numerical Mathematics 54, 171–188
      !! - Mena, H., Ostermann, A., Pfurtscheller, L.-M., Piazzola, C. (2018). "Numerical low-rank 
      !!   approximation of matrix differential equations", Journal of Computational and Applied Mathematics,
      !!   340, 602-614
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_linop_rdp),              intent(inout) :: A
      !! Linear operator
      class(abstract_vector_rdp),             intent(in)    :: B(:)
      !! Low-Rank inhomogeneity.
      real(wp),                               intent(in)    :: Tend
      !! Integration time horizon.
      real(wp),                               intent(inout) :: tau
      !! Desired time step. The avtual time-step will be computed such as to reach Tend in an integer number
      !! of steps.
      integer,                                intent(in)    :: torder
      !! Order of time integration. Only 1st (Lie splitting) and 2nd (Strang splitting) orders are implemented.
      integer,                                intent(out)   :: info
      !! Information flag
      procedure(abstract_exptA_rdp), optional               :: exptA
      !! Routine for computation of the exponential propagator (default: Krylov-based exponential operator).
      logical,                       optional, intent(in)   :: iftrans
      logical                                               :: trans
      !! Determine whether \(\mathbf{A}\) (default `.false.`) or \( \mathbf{A}^T\) (`.true.`) is used.
      logical,                       optional, intent(in)   :: ifverb
      logical                                               :: verbose
      !! Toggle verbosity
      logical,                       optional, intent(in)   :: ifrk
      logical                                               :: rankadaptive
      !! Toggle rank-adaptivity

      ! Internal variables   
      integer                                               :: istep, nsteps, iostep
      real(wp)                                              :: T
      integer                                               :: mode
      procedure(abstract_exptA_rdp), pointer                :: p_exptA => null()
      ! 'timer' to disable rank reduction
      integer                                               :: rk_reduction_lock
      ! global error estimate
      real(wp)                                              :: El
      ! local error estimate
      real(wp)                                              :: err_est
      ! inteval for recomputation of the local error estimate
      integer                                               :: M
      real(wp)                                              :: tol
      character*128                                         :: msg

      ! Optional arguments
      trans        = optval(iftrans, .false.)
      verbose      = optval(ifverb, .false.)
      rankadaptive = optval(ifrk, .false.)
      !tolerance    = optval(tol, 1e-8_wp)
      if (present(exptA)) then
         p_exptA => exptA
      else
         p_exptA => k_exptA_rdp
      end if

      ! Initialize
      T                 = 0.0_wp
      rk_reduction_lock = 10
      M                 = 10

      ! Compute number of steps
      nsteps = floor(Tend/tau)
      if (verbose) write(*,*) 'DLRA integration: nsteps', nsteps

      iostep = nsteps/10
      if ( iostep .eq. 0 ) then
         iostep = 10
      endif

      if ( torder .eq. 1 ) then 
         mode = 1
      else if ( torder .eq. 2 ) then
         mode = 2
      else if ( torder .gt. 2 ) then
         write(msg, *) "Time-integration order for the operator splitting of d > 2 &
                      & requires adjoint solves and is not implemented. Resetting torder = 2." 
         call logger%log_message(trim(msg), module=this_module, procedure='DLRA')
         mode = 2
      else 
         write(msg, *) "Invalid time-integration order specified: ", mode
         call stop_error(trim(msg), module=this_module, &
                           & procedure='numerical_low_rank_splitting_lyapunov_integrator_rdp')
      endif

      ! determine initial rank if rank-adaptive
      if (rankadaptive) then
         call set_initial_rank(X, A, B, tau, mode, p_exptA, trans, 1e-6_wp, verbose=.false.)
         err_est = 0.0_wp
         El      = 0.0_wp
         call compute_splitting_error(err_est, X, A, B, tau, mode, exptA, trans)
         tol = err_est / sqrt(256_wp - real(X%rk + 1))
         if (verbose) then
            write(msg, *) 'Initialization complete: rk = ', X%rk, ', local error estimate: ', tol
            call logger%log_message(trim(msg), module=this_module, procedure='RA-DLRA')
         end if
      end if

      dlra : do istep = 1, nsteps
         ! dynamical low-rank approximation solver
         if (rankadaptive) then
            call DLRA_rank_adaptive_lyapunov_step_rdp(X, A, B, tau, mode, info, rk_reduction_lock, p_exptA, trans, verbose, tol)
            
            if ( mod(istep, M) .eq. 0 ) then
               call compute_splitting_error(err_est, X, A, B, tau, mode, exptA, trans)
               El = El + err_est
               tol = min(El / sqrt(256_wp - real(X%rk + 1)), 1e-5_wp)
               if (verbose) then
                  write(msg, '(3X,I3,A,E8.2)') istep, ': recomputed error estimate: ', tol
                  call logger%log_message(trim(msg), module=this_module, procedure='RA-DLRA')
               end if
            else
               El = El + err_est
            end if
            !
         else
            call numerical_low_rank_splitting_lyapunov_step_rdp(X, A, B, tau, mode, info, p_exptA, trans, verbose)
         end if

         T = T + tau
         ! here we can do some checks such as whether we have reached steady state
         if ( mod(istep,iostep) .eq. 0 ) then
            if (verbose) then
               write(msg, '(3X,I3,A,F6.3)') istep, " steps of DLRA computed. T = ", T
               call logger%log_message(trim(msg), module=this_module, procedure='DLRA-main')
            endif
         endif
      enddo dlra

      if (allocated(U1)) deallocate(U1)
      if (allocated(Uwrk)) deallocate(Uwrk)
      if (allocated(BBTU)) deallocate(BBTU)

      return
   end subroutine numerical_low_rank_splitting_lyapunov_integrator_rdp

   !-----------------------------
   !-----     UTILITIES     -----
   !-----------------------------

   subroutine numerical_low_rank_splitting_lyapunov_step_rdp(X, A, B, tau, mode, info, exptA, trans, verbose)
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_linop_rdp),              intent(inout) :: A
      !! Linear operator
      class(abstract_vector_rdp),             intent(in)    :: B(:)
      !! Low-Rank inhomogeneity.
      real(wp),                               intent(in)    :: tau
      !! Time step.
      integer,                                intent(in)    :: mode
      !! TIme integration mode. Only 1st (Lie splitting - mode 1) and 2nd (Strang splitting - mode 2) orders are implemented.
      integer,                                intent(out)   :: info
      !! Information flag
      procedure(abstract_exptA_rdp)                         :: exptA
      !! Routine for computation of the exponential propagator (default: Krylov-based exponential operator).
      logical,                                intent(in)    :: trans
      !! Determine whether \(\mathbf{A}\) (default `.false.`) or \( \mathbf{A}^T\) (`.true.`) is used.
      logical,                                intent(in)    :: verbose
      !! Verbosity
      
      ! Internal variables
      integer                                               :: istep, nsteps
      character*128                                         :: msg

      select case (mode)
      case (1)
         ! Lie-Trotter splitting
         call M_forward_map(         X, A, tau, info, exptA, trans)
         call G_forward_map_lyapunov(X, B, tau, info)
      case (2) 
         ! Strang splitting
         call M_forward_map(         X, A, 0.5*tau, info, exptA, trans)
         call G_forward_map_lyapunov(X, B,     tau, info)
         call M_forward_map(         X, A, 0.5*tau, info, exptA, trans)
      end select

      return
   end subroutine numerical_low_rank_splitting_lyapunov_step_rdp

   subroutine M_forward_map_rdp(X, A, tau, info, exptA, iftrans)
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_linop_rdp),              intent(inout) :: A
      !! Linear operator.
      real(wp),                               intent(in)    :: tau
      !! Time step.
      integer,                                intent(out)   :: info
      !! Information flag
      procedure(abstract_exptA_rdp)                         :: exptA
      !! Routine for computation of the exponential pabstract_vector),  ropagator (default: Krylov-based exponential operator).
      logical, optional,                      intent(in)    :: iftrans
      logical                                               :: trans
      !! Determine whether \(\mathbf{A}\) (default `.false.`) or \( \mathbf{A}^T\) (`.true.`) is used.

      ! Internal variables
      class(abstract_vector_rdp),             allocatable   :: exptAU    ! scratch basis
      real(wp),                               allocatable   :: R(:,:)    ! QR coefficient matrix
      integer,                                allocatable   :: perm(:)   ! Permutation vector
      real(wp),                               allocatable   :: wrk(:,:)
      integer                                               :: i, rk

      ! Optional argument
      trans = optval(iftrans, .false.)

      rk = X%rk
      allocate(R(rk,rk));   R = 0.0_wp 
      allocate(perm(rk));     perm = 0
      allocate(wrk(rk,rk)); wrk = 0.0_wp

      ! Apply propagator to initial basis
      allocate(exptAU, source=X%U(1)); call exptAU%zero()
      do i = 1, rk
         call exptA(exptAU, A, X%U(i), tau, info, trans)
         call X%U(i)%axpby(0.0_wp, exptAU, 1.0_wp) ! overwrite old solution
      end do
      ! Reorthonormalize in-place
      call qr(X%U(:rk), R, perm, info)
      call check_info(info, 'qr_pivot', module=this_module, procedure='M_forward_map_rdp')
      ! Update low-rank fcators
      call apply_inverse_permutation_matrix(R, perm)
      ! Update coefficient matrix
      wrk = matmul(X%S(:rk,:rk), transpose(R))
      X%S(:rk,:rk) = matmul(R, wrk)

      return
   end subroutine M_forward_map_rdp

   subroutine G_forward_map_lyapunov_rdp(X, B, tau, info)
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_vector_rdp),             intent(in)    :: B(:)
      !! Low-Rank inhomogeneity.
      real(wp),                               intent(in)    :: tau
      !! Time step.
      integer,                                intent(out)   :: info
      !! Information flag.

      ! Internal variables
      class(abstract_vector_rdp),             allocatable   :: U1(:)
      class(abstract_vector_rdp),             allocatable   :: BBTU(:)
      integer                                               :: rk, rkmax

      rk = X%rk
      rkmax = size(X%U)
      if (.not. allocated(U1))   allocate(U1(  rkmax), source=X%U(1))
      if (.not. allocated(BBTU)) allocate(BBTU(rkmax), source=X%U(1))
      call zero_basis(U1); call zero_basis(BBTU)

      call K_step_lyapunov(X, U1(:rk), BBTU(:rk), B, tau, info)
      call S_step_lyapunov(X, U1(:rk), BBTU(:rk),    tau, info)
      call L_step_lyapunov(X, U1(:rk),             B, tau, info)
      
      ! Copy updated low-rank factors to output
      call copy_basis(X%U(:rk), U1(:rk))
               
      return
   end subroutine G_forward_map_lyapunov_rdp

   subroutine K_step_lyapunov_rdp(X, U1, BBTU, B, tau, info)
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_vector_rdp),             intent(out)   :: U1(:)
      !! Intermediate low-rank factor.
      class(abstract_vector_rdp),             intent(out)   :: BBTU(:)
      !! Precomputed application of the inhomogeneity.
      class(abstract_vector_rdp),             intent(in)    :: B(:)
      !! Low-Rank inhomogeneity.
      real(wp),                               intent(in)    :: tau
      !! Time step.
      integer,                                intent(out)   :: info
      !! Information flag.

      ! Internal variables
      integer,                                allocatable   :: perm(:)   ! Permutation vector
      integer                                               :: rk, rkmax

      info = 0

      rk = X%rk
      rkmax = size(X%U)
      if (.not. allocated(Swrk)) allocate(Swrk(rkmax,rkmax)); Swrk = 0.0_wp
      allocate(perm(rk)); perm = 0
      block
         class(abstract_vector_rdp), allocatable :: Xwrk(:)
         call linear_combination(Xwrk, X%U(:rk), X%S(:rk,:rk))             ! K0
         call copy_basis(U1, Xwrk)
      end block
      call apply_outerprod(BBTU, B, X%U(:rk))   ! Kdot
      ! Construct intermediate solution U1
      call axpby_basis(U1, 1.0_wp, BBTU(:rk), tau)   ! K0 + tau*Kdot
      ! Orthonormalize in-place
      call qr(U1(:rk), Swrk(:rk,:rk), perm, info)
      call check_info(info, 'qr_pivot', module=this_module, procedure='K_step_Lyapunov_rdp')
      call apply_inverse_permutation_matrix(Swrk(:rk,:rk), perm)
      X%S(:rk,:rk) = Swrk(:rk,:rk)

      return
   end subroutine K_step_lyapunov_rdp

   subroutine S_step_lyapunov_rdp(X, U1, BBTU, tau, info)
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_vector_rdp),             intent(in)    :: U1(:)
      !! Intermediate low-rank factor.
      class(abstract_vector_rdp),             intent(in)    :: BBTU(:)
      !! Precomputed application of the inhomogeneity.
      real(wp),                               intent(in)    :: tau
      !! Time step.
      integer,                                intent(out)   :: info
      !! Information flag.

      ! Internal variables
      integer                                               :: rk, rkmax

      info = 0

      rk = X%rk
      rkmax = size(X%U)
      if (.not. allocated(Swrk)) allocate(Swrk(rkmax,rkmax)); Swrk = 0.0_wp
      call innerprod(Swrk(:rk,:rk), U1, BBTU)          ! - Sdot
      ! Construct intermediate coefficient matrix
      X%S(:rk,:rk) = X%S(:rk,:rk) - tau*Swrk(:rk,:rk)

      return
   end subroutine S_step_lyapunov_rdp

   subroutine L_step_lyapunov_rdp(X, U1, B, tau, info)
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_vector_rdp),             intent(in)    :: U1(:)
      !! Intermediate low-rank factor (from K step).
      class(abstract_vector_rdp),             intent(in)    :: B(:)
      !! Low-Rank inhomogeneity.
      real(wp),                               intent(in)    :: tau
      !! Time step.
      integer,                                intent(out)   :: info
      !! Information flag.

      ! Internal variables
      integer                                               :: rk, rkmax

      info = 0

      rk = X%rk
      if (.not. allocated(Uwrk)) allocate(Uwrk(rkmax), source=X%U(1)); call zero_basis(Uwrk)

      block
         class(abstract_vector_rdp), allocatable :: Xwrk(:)
         call linear_combination(Xwrk, X%U(:rk), transpose(X%S(:rk,:rk)))  ! L0.T
         call copy_basis(Uwrk(:rk), Xwrk)
      end block
      call apply_outerprod(X%U(:rk), B, U1)       ! Ldot.T
      ! Construct solution L1.T
      call axpby_basis(Uwrk(:rk), 1.0_wp, X%U(:rk), tau)
      ! Update coefficient matrix
      call innerprod(X%S(:rk,:rk), Uwrk(:rk), U1)

      return
   end subroutine L_step_lyapunov_rdp

   !-----------------------------
   !
   !     RANK-ADAPTIVE PSI 
   !
   !-----------------------------

   subroutine DLRA_rank_adaptive_lyapunov_step_rdp(X, A, B, tau, mode, info, rk_reduction_lock, exptA, trans, verbose, tol)
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_linop_rdp),              intent(inout) :: A
      !! Linear operator
      class(abstract_vector_rdp),             intent(in)    :: B(:)
      !! Low-Rank inhomogeneity.
      real(wp),                               intent(in)    :: tau
      !! Time step.
      integer,                                intent(in)    :: mode
      !! TIme integration mode. Only 1st (Lie splitting - mode 1) and 2nd (Strang splitting - mode 2) orders are implemented.
      integer,                                intent(out)   :: info
      !! Information flag
      integer,                                intent(inout) :: rk_reduction_lock
      !! 'timer' to disable rank reduction
      procedure(abstract_exptA_rdp)                         :: exptA
      !! Routine for computation of the exponential propagator (default: Krylov-based exponential operator).
      logical,                                intent(in)    :: trans
      !! Determine whether \(\mathbf{A}\) (default `.false.`) or \( \mathbf{A}^T\) (`.true.`) is used.
      logical,                                intent(in)    :: verbose
      !! Toggle verbosity
      real(wp),                     optional, intent(in)    :: tol
      real(wp)                                              :: tolerance
      !! tolerance for the smallest singular value (this value is only used if ifrk = .true.)
      
      ! Internal variables
      integer                                               :: istep, rk, irk
      logical                                               :: accept_step, found
      real(wp),                               allocatable   :: coef(:)
      real(wp)                                              :: norm
      character*128                                         :: msg

      integer, parameter                                    :: max_step = 5  ! might not be needed

      ! optional arguments
      tolerance = optval(tol, 1e-6_wp)

      ! Allocate memory for SVD
      if (.not.allocated( Usvd)) allocate( Usvd(size(X%U),size(X%U)))
      if (.not.allocated( ssvd)) allocate( ssvd(size(X%U)))
      if (.not.allocated(VTsvd)) allocate(VTsvd(size(X%U),size(X%U)))

      ! ensure that we are integrating one more rank than we use for approximation
      X%rk = X%rk + 1
      rk = X%rk ! this is only to make the code more readable
      
      accept_step = .false.
      istep = 1
      do while ( .not. accept_step .and. istep < max_step )
         ! run a regular step
         call numerical_low_rank_splitting_lyapunov_step_rdp(X, A, B, tau, mode, info, exptA, trans, verbose)
         ! compute singular values of X%S
         Usvd = 0.0_wp; ssvd = 0.0_wp; VTsvd = 0.0_wp
         call svd(X%S(:rk,:rk), Usvd(:rk,:rk), ssvd(:rk), VTsvd(:rk,:rk))
         found = .false.
         tol_chk: do irk = 1, rk
            if ( ssvd(irk) < tolerance ) then
               found = .true.
               exit tol_chk
            end if
         end do tol_chk
         if (.not. found) irk = irk - 1
         
         ! choose action
         if (.not. found) then ! none of the singular values is below tolerance
            ! increase rank and run another step
            if (rk == size(X%U)) then ! cannot increase rank without reallocating X%U and X%S
               write(msg, *) 'Cannot increase rank, rkmax is reached. Increase rkmax and restart!'
               call stop_error(trim(msg), module=this_module, procedure='DLRA_rank_adaptive_lyapunov_step_rdp')
            else
               if (verbose) then
                  write(msg,'(6X,A,I3)') '--> increase rank to ', rk + 1
                  call logger%log_message(trim(msg), module=this_module, procedure='RA-DLRA')
               end if
               X%rk = X%rk + 1
               rk = X%rk ! this is only to make the code more readable
               ! set coefficients to zero (for redundancy)
               X%S(:rk, rk) = 0.0_wp 
               X%S( rk,:rk) = 0.0_wp
               ! add random vector ...
               call X%U(rk)%rand(.false.)
               ! ... and orthonormalize
               call orthogonalize_against_basis(X%U(rk), X%U(:rk-1), info, if_chk_orthonormal=.false.)
               call check_info(info, 'orthogonalize_against_basis', module=this_module, &
                                 & procedure='DLRA_rank_adaptive_lyapunov_step_rdp')
               call X%U(rk)%scal(1.0_wp / X%U(rk)%norm())

               rk_reduction_lock = 10 ! avoid rank oscillations

            end if
         else ! the rank of the solution is sufficient
            accept_step = .true.

            if (irk /= rk .and. rk_reduction_lock == 0) then ! we should decrease the rank
               ! decrease rank
               
               ! rotate basis onto principal axes
               block
                  class(abstract_vector_rdp), allocatable :: Xwrk(:)
                  call linear_combination(Xwrk, X%U(:rk), Usvd(:rk,:rk))
                  call copy_basis(X%U(:rk), Xwrk)
               end block
               X%S(:rk,:rk) = diag(ssvd(:rk))

               rk = max(irk, rk - 2)  ! reduce by at most 2

               if (verbose) then
                  write(msg, '(6X,A,I3)') '--> decrease rank to ', rk
                  call logger%log_message(trim(msg), module=this_module, procedure='RA-DLRA')
               end if
            end if
            
         end if ! found
         istep = istep + 1
      end do ! while .not. accept_step
      if (verbose) then
         write(msg,'(A,I3,A,I2,A,E14.8,A,I2)') 'rk = ', X%rk-1, ':     s_', irk,' = ', &
                                                   & ssvd(irk), ', rank_lock: ', rk_reduction_lock
         call logger%log_message(trim(msg), module=this_module, procedure='RA-DLRA')
      end if

      ! decrease rk_reduction_lock
      if (rk_reduction_lock > 0) rk_reduction_lock = rk_reduction_lock - 1
      
      ! reset to the rank of the approximation which we use outside of the integrator
      X%rk = rk - 1      

      return
   end subroutine DLRA_rank_adaptive_lyapunov_step_rdp

   subroutine set_initial_rank(X, A, B, tau, mode, exptA, trans, tol, verbose, rk_init, nsteps)
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_linop_rdp),              intent(inout) :: A
      !! Linear operator
      class(abstract_vector_rdp),             intent(in)    :: B(:)
      !! Low-Rank inhomogeneity.
      real(wp),                               intent(in)    :: tau
      !! Time step.
      integer,                                intent(in)    :: mode
      !! TIme integration mode. Only 1st (Lie splitting - mode 1) and 2nd (Strang splitting - mode 2) orders are implemented.
      procedure(abstract_exptA_rdp)                         :: exptA
      !! Routine for computation of the exponential propagator (default: Krylov-based exponential operator).
      logical,                                intent(in)    :: trans
      !! Determine whether \(\mathbf{A}\) (default `.false.`) or \( \mathbf{A}^T\) (`.true.`) is used.
      real(wp),                               intent(in)    :: tol
      !! Tolerance on the last singular value to determine rank
      logical,                                intent(in)    :: verbose
      !! verbosity
      integer,                      optional, intent(in)    :: rk_init
      !! Smallest tested rank
      integer,                      optional, intent(in)    :: nsteps
      integer                                               :: n
      !! Number of steps to run before checking the singular values

      ! internal
      integer                                               :: i, irk, info, rkmax
      class(abstract_vector_rdp),               allocatable :: Utmp(:)
      real(wp),                                 allocatable :: Stmp(:,:)
      logical                                               :: found, accept_rank
      character*128                                         :: msg

      ! optional arguments
      X%rk = optval(rk_init, 1)
      n = optval(nsteps, 5)
      rkmax = size(X%U)

      ! Allocate memory for SVD
      if (.not.allocated( Usvd)) allocate( Usvd(rkmax,rkmax))
      if (.not.allocated( ssvd)) allocate( ssvd(rkmax))
      if (.not.allocated(VTsvd)) allocate(VTsvd(rkmax,rkmax))

      info = 0
      accept_rank = .false.

      ! save initial condition
      allocate(Utmp(rkmax), source=X%U)
      allocate(Stmp(rkmax,rkmax)); Stmp = X%S
      
      do while (.not. accept_rank .and. X%rk <= rkmax)
         ! run integrator
         do i = 1,n
            call numerical_low_rank_splitting_lyapunov_step_rdp(X, A, B, tau, mode, info, exptA, trans, verbose)
         end do

         ! check if singular values are resolved
         Usvd = 0.0_wp; ssvd = 0.0_wp; VTsvd = 0.0_wp
         call svd(X%S(:X%rk,:X%rk), Usvd(:X%rk,:X%rk), ssvd(:X%rk), VTsvd(:X%rk,:X%rk))
         found = .false.
         tol_chk: do irk = 1, X%rk
            if ( ssvd(irk) < tol ) then
               found = .true.
               exit tol_chk
            end if
         end do tol_chk
         if (.not. found) irk = irk - 1
         if (verbose) then
            write(msg,'(A,I2,A,E8.2)') '   rk = ', X%rk, ' s_r =', ssvd(X%rk)
            call logger%log_message(trim(msg), module=this_module, procedure='numerical_low_rank_splitting_lyapunov_integrator_rdp')
         end if
         if (found) then
            accept_rank = .true.
            X%rk = irk
            write(msg,'(A,I2,A,E10.4)') '   Accpeted rank: r = ', X%rk-1, ',     s_{r+1} = ', ssvd(X%rk)
            call logger%log_message(trim(msg), module=this_module, procedure='set_initial_rank')
         else
            X%rk = 2*X%rk
         end if
         
         ! reset initial conditions
         call copy_basis(X%U, Utmp)
         X%S = Stmp
      end do

      if (X%rk > rkmax) then
         write(msg, *) 'Maximum rank reached but singular values are not converged. Increase rkmax and restart.'
         call stop_error(trim(msg), module=this_module, procedure='set_initial_rank')
      end if

      ! reset to the rank of the approximation which we use outside of the integrator 
      X%rk = X%rk - 1

   end subroutine set_initial_rank

   subroutine compute_splitting_error(err_est, X, A, B, tau, mode, exptA, trans)
      real(wp),                               intent(out)   :: err_est
      !! Estimation of the splitting error
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_linop_rdp),              intent(inout) :: A
      !! Linear operator
      class(abstract_vector_rdp),             intent(in)    :: B(:)
      !! Low-Rank inhomogeneity.
      real(wp),                               intent(in)    :: tau
      !! Time step.
      integer,                                intent(in)    :: mode
      !! TIme integration mode. Only 1st (Lie splitting - mode 1) and 2nd (Strang splitting - mode 2) orders are implemented.
      procedure(abstract_exptA_rdp)                         :: exptA
      !! Routine for computation of the exponential propagator (default: Krylov-based exponential operator).
      logical,                                intent(in)    :: trans
      !! Determine whether \(\mathbf{A}\) (default `.false.`) or \( \mathbf{A}^T\) (`.true.`) is used.
      
      ! internals
      class(abstract_vector_rdp),               allocatable :: Utmp(:)
      real(wp),                                 allocatable :: Stmp(:,:)
      class(abstract_vector_rdp),               allocatable :: U1(:), Vperp(:)
      real(wp),                                 allocatable :: S1(:,:)
      real(wp),                                 allocatable :: Vt1(:,:), Vt2(:,:)
      real(wp),                                 allocatable :: Sdiff(:,:)
      integer                                               :: info

      ! svd
      real(wp),                                 allocatable :: ssvd(:)
      real(wp),                                 allocatable :: Usvd(:,:), VTsvd(:,:)

      ! scratch
      allocate(Vt1(X%rk,X%rk)); Vt1 = 0.0_wp
      allocate(Vt2(X%rk,X%rk)); Vt2 = 0.0_wp
      allocate(Sdiff(2*X%rk,2*X%rk)); Sdiff = 0.0_wp

      ! save curret state
      allocate(Utmp(X%rk), source=X%U(:X%rk))
      allocate(Stmp(X%rk,X%rk)); Stmp = X%S(:X%rk,:X%rk)

      ! tau step
      call numerical_low_rank_splitting_lyapunov_step_rdp(X, A, B, tau, mode, info, exptA, trans, verbose=.false.)
      ! save result
      allocate(U1(X%rk), source=X%U(:X%rk))
      allocate(S1(X%rk,X%rk)); S1 = X%S(:X%rk,:X%rk)

      ! reset curret state
      call copy_basis(X%U(:X%rk), Utmp)
      X%S(:X%rk,:X%rk) = Stmp
      ! tau/2 steps
      call numerical_low_rank_splitting_lyapunov_step_rdp(X, A, B, 0.5*tau, mode, info, exptA, trans, verbose=.false.)
      call numerical_low_rank_splitting_lyapunov_step_rdp(X, A, B, 0.5*tau, mode, info, exptA, trans, verbose=.false.)

      ! compute common basis
      allocate(Vperp(X%rk), source=X%U(:X%rk))  ! V_perp = V
      call orthogonalize_against_basis(Vperp, U1, info, if_chk_orthonormal=.false., beta=Vt1) ! beta = Vt1.T
      call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='compute_splitting_error')
      call qr(Vperp, Swrk, info)
      call check_info(info, 'qr', module=this_module, procedure='compute_splitting_error')
      call innerprod(Vt2, X%U(1:X%rk), Vperp)

      ! project X_2 onto extended basis and construct difference matrix
      Sdiff(      :  X%rk,      :  X%rk) = S1 - matmul(          Vt1,  matmul(X%S(:X%rk,:X%rk), transpose(Vt1)))
      Sdiff(X%rk+1:2*X%rk,      :  X%rk) =    - matmul(transpose(Vt2), matmul(X%S(:X%rk,:X%rk), transpose(Vt1)))
      Sdiff(      :  X%rk,X%rk+1:2*X%rk) =    - matmul(          Vt1,  matmul(X%S(:X%rk,:X%rk),           Vt2))
      Sdiff(X%rk+1:2*X%rk,X%rk+1:2*X%rk) =    - matmul(transpose(Vt2), matmul(X%S(:X%rk,:X%rk),           Vt2))
      
      ! svd
      allocate( Usvd(2*X%rk,2*X%rk))
      allocate( ssvd(2*X%rk))
      allocate(VTsvd(2*X%rk,2*X%rk))
      call svd(Sdiff, Usvd, ssvd, VTsvd)

      ! compute local error based on frobenius norm of difference
      err_est = 2**mode / (2**mode - 1) * sqrt( sum( ssvd ** 2 ) )

      ! reset curret state
      call copy_basis(X%U(:X%rk), Utmp)
      X%S(:X%rk,:X%rk) = Stmp

   end subroutine compute_splitting_error

end module LightROM_LyapunovSolvers
