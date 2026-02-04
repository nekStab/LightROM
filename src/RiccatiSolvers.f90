module LightROM_RiccatiSolvers
   !! This module provides the implementation of the Krylov-based solvers for the Differential Riccati
   !! equation based on the dynamic low-rank approximation and operator splitting.
   ! Standard Library
   use stdlib_optval, only : optval
   use stdlib_linalg, only : eye
   ! LightKrylov modules
   use LightKrylov
   use LightKrylov, only: wp => dp
   use LightKrylov_Logger
   use LightKrylov_AbstractVectors
   use LightKrylov_ExpmLib
   use Lightkrylov_BaseKrylov
   ! LightROM modules
   use LightROM_Utils
   use LightROM_LoggerUtils
   use LightROM_Timing, only: lr_timer => global_lightROM_timer, time_lightROM
   use LightROM_LyapunovUtils
   use LightROM_LyapunovSolvers, only : M_forward_map
   use LightROM_RiccatiUtils
   use LightROM_AbstractLTIsystems

   implicit none

   ! global scratch arrays
   class(abstract_vector_rdp),  allocatable   :: Uwrk0(:)
   class(abstract_vector_rdp),  allocatable   :: Uwrk1(:)
   class(abstract_vector_rdp),  allocatable   :: U1(:)
   class(abstract_vector_rdp),  allocatable   :: QU(:)
   real(wp),                    allocatable   :: Swrk0(:,:)
   real(wp),                    allocatable   :: Swrk1(:,:)
   ! global scratch arrays for the predictor step
   class(abstract_vector_rdp),  allocatable   :: U0(:)
   class(abstract_vector_rdp),  allocatable   :: T0(:)
   class(abstract_vector_rdp),  allocatable   :: Ut(:)
   class(abstract_vector_rdp),  allocatable   :: Tt(:)
   real(wp),                    allocatable   :: S0(:,:)
   ! global scratch arrays SVD
   real(wp),                    allocatable   :: ssvd(:)
   real(wp),                    allocatable   :: Usvd(:,:), VTsvd(:,:)

   private 
   ! module name
   private :: this_module
   character(len=*), parameter :: this_module = 'LR_RiccSolvers'
   character(len=*), parameter :: logfile_basename = 'Riccati_SVD'
   integer :: RiccSolver_counter = 0

   public :: projector_splitting_DLRA_riccati_integrator
   public :: G_forward_map_riccati
   public :: K_step_riccati
   public :: S_step_riccati
   public :: L_step_riccati
   public :: reset_riccati_solver

   interface projector_splitting_DLRA_riccati_integrator
      module procedure projector_splitting_DLRA_riccati_integrator_rdp
   end interface

   interface G_forward_map_riccati
      module procedure G_forward_map_riccati_rdp
   end interface

   interface K_step_riccati
      module procedure K_step_riccati_rdp
   end interface

   interface S_step_riccati
      module procedure S_step_riccati_rdp
   end interface

   interface L_step_riccati
      module procedure L_step_riccati_rdp
   end interface

   contains

   subroutine projector_splitting_DLRA_riccati_integrator_rdp(X, A, B, CT, Qc, Rinv, Tend, tau, info, &
                                                               & exptA, iftrans, options)
      !! Main driver for the numerical integrator for the matrix-valued differential Riccati equation of the form
      !!
      !!    $$\dot{\mathbf{X}} = \mathbf{A} \mathbf{X} + \mathbf{X} \mathbf{A}^T + \mathbf{C}^T \mathbf{Q} \mathbf{C} - \mathbf{X} \mathbf{B} \mathbf{R}^{-1} \mathbf{B}^T \mathbf{X} $$
      !!
      !! where \( \mathbf{A} \) is a (n x n) Hurwitz matrix, \( \mathbf{X} \) is SPD and 
      !! \( \mathbf{B} \) and \( \mathbf{C}^T \) are low-rank matrices.
      !!
      !! Since \( \mathbf{A} \) is Hurwitz, the equations converges to steady state for  \( t \to \infty \), 
      !! which corresponds to the associated algebrait Riccati equation of the form
      !!
      !!    $$\mathbf{0} = \mathbf{A} \mathbf{X} + \mathbf{X} \mathbf{A}^T + \mathbf{C}^T \mathbf{Q} \mathbf{C} - \mathbf{X} \mathbf{B} \mathbf{R}^{-1} \mathbf{B}^T \mathbf{X} $$
      !!
      !! The algorithm is based on four main ideas:
      !!
      !! - Dynamic Low-Rank Approximation (DLRA). DLRA is a method for the solution of general matrix differential 
      !!   equations proposed by Nonnenmacher & Lubich (2007) which seeks to integrate only the leading low-rank 
      !!   factors of the solution to a large system by updating an appropriate matrix factorization. The time-integration
      !!   is achieved by splitting the step into three sequential substeps, each updating a part of the factorization
      !!   taking advantage of and maintaining the orthogonality of the left and right low-rank bases of the factorization.
      !! - Projector-Splitting Integration (PSI). The projector-splitting scheme proposed by Lubich & Oseledets (2014) 
      !!   for the solution of DLRA splits the right-hand side of the differential equation into a linear stiff part 
      !!   that is integrated exactly and a (possibly non-linear) non-stiff part which is integrated numerically. 
      !!   The two operators are then composed to obtain the integrator for the full differential equation.
      !!   The advantage of the projector splitting integration is that it maintains orthonormality of the basis
      !!   of the low-rank approximation to the solution without requiring SVDs of the full matrix.                                                                     
      !! - The third element is the application of the general framework of projector-splitting integration for 
      !!   dynamical low-rank approximation to the Riccati equations by Mena et al. (2018). As the solutions
      !!   to the Riccati equation are by construction SPD, this fact can be taken advantage of to reduce the 
      !!   computational cost of the integration and, in particular, doing away with one QR factorization per timestep
      !!   while maintaining symmetry of the resulting matrix factorization.
      !! - The final element is the addition of the capability of dyanmic rank adaptivity for the projector-splitting
      !!   integrator proposed by Hochbruck et al. (2023). At the cost of integrating a supplementary solution vector, 
      !!   the rank of the solution is dynamically adapted to ensure that the corresponding additional singular value
      !!   stays below a chosen threshold.
      !!
      !! **Algorithmic Features**
      !! 
      !! - Separate integration of the stiff inhomogeneous part of the Riccati equation and the non-stiff inhomogeneity
      !! - Rank preserving time-integration that maintains orthonormality of the factorization basis
      !! - Alternatively, dynamical rank-adaptivity based on the instantaneous singular values
      !! - The stiff part of the problem is solved using a time-stepper approach to approximate 
      !!   the action of the exponential propagator
      !!
      !! **Advantages**
      !!
      !! - Rank of the approximate solution is user defined or chosen adaptively based on the solution
      !! - The integrator is adjoint-free
      !! - The operator of the homogeneous part and the inhomogeneity are not needed explicitly i.e. the algorithm 
      !! is amenable to solution using Krylov methods (in particular for the solution of the stiff part of the problem)
      !! - No SVDs of the full solution are required for this algorithm
      !! - Lie and Strang splitting implemented allowing for first and second order integration in time
      !!
      !! ** Limitations**
      !!
      !! - Rank of the approximate solution is user defined. The appropriateness of this approximation is not considered.
      !!   This does not apply to the rank-adaptive version of the integrator.
      !! - The current implementation does not require an adjoint integrator. This means that the temporal order of the 
      !!   basic operator splitting scheme is limited to 1 (Lie-Trotter splitting) or at most 2 (Strang splitting). 
      !!   Higher order integrators are possible, but require at least some backward integration (via the adjoint) 
      !!   in BOTH parts of the splitting (see Sheng-Suzuki and Goldman-Kaper theorems).
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
      !! - Hochbruck, M., Neher, M., Schrammer, S. (2023). "Rank-adaptive dynamical low-rank integrators for
      !!   first-order and second-order matrix differential equations", BIT Numerical Mathematics 63:9
      class(abstract_sym_low_rank_state_rdp),  intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_linop_rdp),               intent(inout) :: A
      !! Linear operator.
      class(abstract_vector_rdp),              intent(in)    :: B(:)
      !! System input.
      class(abstract_vector_rdp),              intent(in)    :: CT(:)
      !! System output.
      real(wp),                                intent(in)    :: Qc(:,:)
      !! Measurement weights.
      real(wp),                                intent(in)    :: Rinv(:,:)
      !! Inverse of the actuation weights.
      real(wp),                                intent(in)    :: Tend
      !! Integration time horizon. 
      real(wp),                                intent(in)    :: tau
      !! Desired time step. The avtual time-step will be computed such as to reach Tend in an integer number
      !! of steps.
      integer,                                 intent(out)   :: info
      !! Information flag.
      procedure(abstract_exptA_rdp)                          :: exptA
      !! Routine for computation of the exponential propagator (default: Krylov-based exponential operator).
      logical,                       optional, intent(in)    :: iftrans
      logical                                                :: trans
      !! Determine whether \(\mathbf{A}\) (default `.false.`) or \( \mathbf{A}^T\) (`.true.`) is used.
      type(dlra_opts),               optional, intent(in)    :: options
      type(dlra_opts)                                        :: opts
      !! Options for solver configuration

      ! Internal variables
      integer                             :: i, irk
      integer                             :: istep, nsteps, rkmax, chkstep
      integer                             :: rk_reduction_lock   ! 'timer' to disable rank reduction
      real(wp)                            :: tol                 ! current tolerance
      logical                             :: if_lastep
      real(wp), dimension(:), allocatable :: svals, svals_lag
      character(len=128)                  :: msg

      if (time_lightROM()) call lr_timer%start('DLRA_riccati_integrator_rdp')
      RiccSolver_counter = RiccSolver_counter + 1

      ! Optional arguments
      trans = optval(iftrans, .false.)

      ! Options
      if (present(options)) then
         opts = options
      else ! default
         opts = dlra_opts()
      end if

      ! Set tolerance
      tol = opts%tol

      ! Check compatibility of options and determine chk/IO step
      call check_options(chkstep, tau, X, opts)

      ! Initialize
      rk_reduction_lock = 10
      X%is_converged = .false.
      X%time = 0.0_wp
      X%step = 0
      ! Compute number of steps
      if_lastep = .false.
      nsteps = nint(Tend/tau)

      rkmax = size(X%U)
      ! Allocate memory for SVD & lagged fields
      allocate(Usvd(rkmax,rkmax), ssvd(rkmax), VTsvd(rkmax,rkmax))

      call log_message('Initializing Riccati solver', this_module, 'DLRA_main')
      write(msg,'(A,I0,A,F10.8)') 'Integration over ', nsteps, ' steps with dt= ', tau
      call log_information(msg, this_module, 'DLRA_main')
      ! Prepare logfile
      call write_logfile_headers(logfile_basename)

      if ( opts%mode > 2 ) then
         write(msg,'(A)') "Time-integration order for the operator splitting of d > 2 &
                      & requires adjoint solves and is not implemented. Resetting torder = 2." 
         call log_message(msg, this_module, 'DLRA_main')
      else if ( opts%mode < 1 ) then
         write(msg,'(A,I0)') "Invalid time-integration order specified: ", opts%mode
         call stop_error(msg, this_module, 'DLRA_main')
      else if ( opts%mode == 2 ) then
         write(msg,'(A)') "Second order time-integration order currently not implemented."
         call stop_error(msg, this_module, 'DLRA_main')
      endif

      ! determine initial rank if rank-adaptive
      if (opts%if_rank_adaptive) then
         if (.not. X%rank_is_initialised) then
            call log_message('Determine initial rank:', this_module, 'DLRA_main')
            call set_initial_rank_riccati(X, A, B, CT, Qc, Rinv, tau, opts%mode, exptA, trans, tol)
         end if
      end if

      call log_settings(X, Tend, tau, nsteps, opts)
      call log_message('Starting DLRA integration', this_module, 'DLRA_main')

      dlra : do istep = 1, nsteps

         call log_step(X, istep, nsteps)

         ! save lag data defore the timestep
         if ( X%rk > 0 .and. (mod(istep, chkstep) == 0 .or. istep == nsteps) ) then
            svals_lag = svdvals(X%S(:X%rk,:X%rk))
         end if

         ! dynamical low-rank approximation solver
         if (opts%if_rank_adaptive) then
            call rank_adaptive_PS_DLRA_riccati_step_rdp(X, A, B, CT, Qc, Rinv, tau, opts%mode, info, rk_reduction_lock, & 
                                                         & exptA, trans, tol)
         else
            call projector_splitting_DLRA_riccati_step_rdp(X, A, B, CT, Qc, Rinv, tau, opts%mode, info, exptA, trans)
         end if

         ! update time & step counters
         X%time = X%time + tau
         X%step = X%step + 1
         X%tot_time = X%tot_time + tau
         X%tot_step = X%tot_step + 1

         ! here we can do some checks such as whether we have reached steady state
         if ( X%rk > 0 .and. (mod(istep, chkstep) == 0 .or. istep == nsteps) ) then
            svals = svdvals(X%S(:X%rk,:X%rk))
            if (.not. allocated(svals_lag)) allocate(svals_lag(X%rk), source=zero_rdp)
            call log_svals(logfile_basename, X, tau, svals, svals_lag, RiccSolver_counter, istep, nsteps)
            ! Check convergence
            if (istep == nsteps) if_lastep = .true.
            irk = min(size(svals), size(svals_lag))
            X%is_converged = is_converged(X, svals(:irk), svals_lag(:irk), opts, if_lastep)
            if (X%is_converged) then
               write(msg,'(A,I0,A)') "Step ", istep, ": Solution converged!"
               call log_information(msg, this_module, 'DLRA_main')
               exit dlra
            else ! if final step
               if (if_lastep) then
                  write(msg,'(A,I0,A)') "Step ", istep, ": Solution not converged!"
                  call log_information(msg, this_module, 'DLRA_main')
               end if
            end if
         endif
      enddo dlra
      call log_message('Exiting Riccati solver', this_module, 'DLRA_main')
      ! Clean up scratch space
      if (allocated(Uwrk0)) deallocate(Uwrk0,Uwrk1,U1,QU,Swrk0,Swrk1)
      if (allocated(Usvd)) deallocate(Usvd, ssvd, VTsvd)
      if (time_lightROM()) call lr_timer%stop('DLRA_riccati_integrator_rdp')
   end subroutine projector_splitting_DLRA_riccati_integrator_rdp

   !-----------------------
   !-----     PSI     -----
   !-----------------------

   subroutine projector_splitting_DLRA_riccati_step_rdp(X, A, B, CT, Qc, Rinv, tau, mode, info, exptA, trans)
      !! Driver for the time-stepper defining the splitting logic for each step of the the 
      !! projector-splitting integrator
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_linop_rdp),              intent(inout) :: A
      !! Linear operator.
      class(abstract_vector_rdp),             intent(in)    :: B(:)
      !! System input.
      class(abstract_vector_rdp),             intent(in)    :: CT(:)
      !! System output.
      real(wp),                               intent(in)    :: Qc(:,:)
      !! Measurement weights.
      real(wp),                               intent(in)    :: Rinv(:,:)
      !! Inverse of the actuation weights.
      real(wp),                               intent(in)    :: tau
      !! Time step.
      integer,                                intent(in)    :: mode
      !! Order of time integration. Only 1st (Lie splitting) and 2nd (Strang splitting) 
      !! orders are implemented.
      integer,                                intent(out)   :: info
      !! Information flag.
      procedure(abstract_exptA_rdp), optional               :: exptA
      !! Routine for computation of the exponential propagator (default: Krylov-based exponential operator).
      logical,                   optional,    intent(in)    :: trans
      !! Determine whether \(\mathbf{A}\) (default `.false.`) or \( \mathbf{A}^T\) (`.true.`) is used.

      ! Internal variables
      integer                                               :: istep, nsteps, rk

      if (time_lightROM()) call lr_timer%start('DLRA_riccati_step_rdp')

      select case (mode)
      case (1) 
         ! Lie-Trotter splitting
         call M_forward_map        (X, A,                      tau, info, exptA, trans)
         call G_forward_map_riccati(X,    B, CT, Qc, Rinv,     tau, info)
      case (2) 
         ! Strang splitting
         ! Prepare arrays for predictor step
         rk = size(X%U)
         ! precomputation arrays
         if (.not. allocated(U0)) allocate(U0(1:rk), source=X%U(1:rk))
         if (.not. allocated(T0)) allocate(T0(1:rk), source=X%U(1:rk))
         if (.not. allocated(Ut)) allocate(Ut(1:rk), source=X%U(1:rk))
         if (.not. allocated(Tt)) allocate(Tt(1:rk), source=X%U(1:rk))
         if (.not. allocated(S0)) allocate(S0(1:rk,1:rk))
         call zero_basis(U0); call zero_basis(T0); call zero_basis(Ut); call zero_basis(Tt); S0 = 0.0_wp
         ! scratch arrays
         if (.not. allocated(Uwrk0)) allocate(Uwrk0(1:rk), source=X%U(1:rk))
         if (.not. allocated(Swrk0)) allocate(Swrk0(1:rk,1:rk))
         call zero_basis(Uwrk0); Swrk0 = 0.0_wp
         ! second order step
         call M_forward_map        (X, A,                  0.5*tau, info, exptA, trans)
         ! Save current state
         call copy(U0, X%U); S0 = X%S               ! --> save  
         ! Precompute T0
         block
            class(abstract_vector_rdp), allocatable :: Xwrk(:)
            call linear_combination(Xwrk, U0, S0); call copy(Uwrk0, Xwrk) ! K0 = U0 @ S0
            call apply_premult_outerprod_w(Swrk0, U0, Uwrk0, B, Rinv)   ! (U0.T) @ B @ R^(-1) @ B.T @ K0
            call linear_combination(Xwrk, Uwrk0, Swrk0); call copy(T0, Xwrk) ! K0 @ Swrk0
         end block
         ! First order integration
         call G_forward_map_riccati(X,    B, CT, Qc, Rinv,     tau, info, ifpred=.true., T0=T0)
         ! Precompute Tt
         call copy(Ut, X%U)                                               ! --> save
         block
            class(abstract_vector_rdp), allocatable :: Xwrk(:)
            call linear_combination(Xwrk, X%U, X%S); call copy(Uwrk0, Xwrk) ! Kt = Ut @ St
            call apply_premult_outerprod_w(Swrk0, X%U, Uwrk0, B, Rinv)  ! (Ut.T) @ B @ R^(-1) @ B.T @ Kt
            call linear_combination(Xwrk, Uwrk0, Swrk0); call copy(Tt, Xwrk) ! Kt @ Swrk0
         end block
         ! Reset state
         call copy(X%U, U0); X%S = S0
         ! Second order integration
         call G_forward_map_riccati(X,    B, CT, Qc, Rinv,     tau, info, ifpred=.false., &
                                  & T0=T0, Tt=Tt, U0=U0, Ut=Ut)
         call M_forward_map        (X, A,                  0.5*tau, info, exptA, trans)
      end select

      if (time_lightROM()) call lr_timer%stop('DLRA_riccati_step_rdp')
   end subroutine projector_splitting_DLRA_riccati_step_rdp

   !-----------------------------
   !
   !     RANK-ADAPTIVE PSI 
   !
   !-----------------------------

   subroutine rank_adaptive_PS_DLRA_riccati_step_rdp(X, A, B, CT, Qc, Rinv, tau, mode, info, rk_reduction_lock, exptA, trans, tol)
      !! Wrapper for projector_splitting_DLRA_riccati_step_rdp adding the logic for rank-adaptivity
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_linop_rdp),              intent(inout) :: A
      !! Linear operator
      class(abstract_vector_rdp),             intent(in)    :: B(:)
      !! System input.
      class(abstract_vector_rdp),             intent(in)    :: CT(:)
      !! System output.
      real(wp),                               intent(in)    :: Qc(:,:)
      !! Measurement weights.
      real(wp),                               intent(in)    :: Rinv(:,:)
      !! Inverse of the actuation weights.
      real(wp),                               intent(in)    :: tau
      !! Time step.
      integer,                                intent(in)    :: mode
      !! Time integration mode. Only 1st (Lie splitting - mode 1) and 2nd (Strang splitting - mode 2) 
      !! orders are implemented.
      integer,                                intent(out)   :: info
      !! Information flag
      integer,                                intent(inout) :: rk_reduction_lock
      !! 'timer' to disable rank reduction
      procedure(abstract_exptA_rdp)                         :: exptA
      !! Routine for computation of the exponential propagator (default: Krylov-based exponential operator).
      logical,                                intent(in)    :: trans
      !! Determine whether \(\mathbf{A}\) (default `.false.`) or \( \mathbf{A}^T\) (`.true.`) is used.
      real(wp),                               intent(in)    :: tol
      
      ! Internal variables
      integer                                               :: istep, rk, irk, rkmax, ndigits
      logical                                               :: accept_step, found
      real(wp),                               allocatable   :: coef(:)
      real(wp)                                              :: norm
      character(len=256)                                    :: msg, fmt

      integer, parameter                                    :: max_step = 40  ! might not be needed

      ! ensure that we are integrating one more rank than we use for approximation
      X%rk = X%rk + 1
      rk = X%rk ! this is only to make the code more readable
      rkmax = size(X%U)
      ndigits = max(1,ceiling(log10(real(rkmax))))
      
      accept_step = .false.
      istep = 1
      do while ( .not. accept_step .and. istep < max_step )
         ! run a regular step
         call projector_splitting_DLRA_riccati_step_rdp(X, A, B, CT, Qc, Rinv, tau, mode, info, exptA, trans)
         ! compute singular values of X%S
         Usvd = 0.0_wp; ssvd = 0.0_wp; VTsvd = 0.0_wp
         call svd(X%S(:rk,:rk), ssvd(:rk), Usvd(:rk,:rk), VTsvd(:rk,:rk))
         found = .false.
         tol_chk: do irk = 1, rk
            if ( ssvd(irk) < tol ) then
               found = .true.
               exit tol_chk
            end if
         end do tol_chk
         !if (.not. found) irk = irk - 1
         
         ! choose action
         if (.not. found) then ! none of the singular values is below tolerance
            ! increase rank and run another step
            if (rk == rkmax) then ! cannot increase rank without reallocating X%U and X%S
               write(msg,'(A,I0,A,A)') 'Cannot increase rank, rkmax = ', rkmax, ' is reached. ', &
                        & 'Increase rkmax and restart!'
               call stop_error(msg, this_module, 'rank_adaptive_PS_DLRA_lyapunov_step_rdp')
            else
               write(fmt,'("(A,I3,A,I",I0,".",I0,",A,E14.8)")') ndigits, ndigits
               write(msg,fmt) 'rk= ', rk, ', s_', rk,' = ', ssvd(rk)
               call log_information(msg, this_module, 'DLRA_main')
               write(msg,'(A,I0)') 'Rank increased to rk= ', rk + 1
               call log_message(msg, this_module, 'DLRA_main')
               
               X%rk = X%rk + 1
               rk = X%rk ! this is only to make the code more readable
               ! set coefficients to zero (for redundancy)
               X%S(:rk, rk) = 0.0_wp 
               X%S( rk,:rk) = 0.0_wp
               ! add random vector ...
               call X%U(rk)%rand(.false.)
               ! ... and orthonormalize
               call orthogonalize_against_basis(X%U(rk), X%U(:rk-1), info, if_chk_orthonormal=.false.)
               call check_info(info, 'orthogonalize_against_basis', this_module, &
                                 & 'rank_adaptive_PS_DLRA_lyapunov_step_rdp')
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
                  call copy(X%U(:rk), Xwrk)
               end block
               X%S(:rk,:rk) = diag(ssvd(:rk))

               rk = max(irk, rk - 2)  ! reduce by at most 2

               write(msg, '(A,I0)') 'Rank decreased to rk= ', rk
               call log_message(msg, this_module, 'DLRA_main')
            end if
            
         end if ! found
         istep = istep + 1
      end do ! while .not. accept_step

      if (.not. accept_step .and. istep == max_step) then
         write(msg,'(A,I0,A,2(A,E11.4))') 'Rank increased ', max_step, ' times in a single step without ', &
               & 'reaching the desired tolerance on the singular values. s_{k+1} = ', ssvd(irk), ' > ', tol
         call stop_error(msg, this_module, 'DLRA_main')
      end if

      write(fmt,'("(A,I3,A,I",I0,".",I0,",A,E14.8,A,I2)")') ndigits, ndigits
      write(msg,fmt) 'rk = ', rk-1, ':     s_', irk,' = ', &
               & ssvd(irk), ',     lock: ', rk_reduction_lock
      call log_information(msg, this_module, 'DLRA_main')

      ! decrease rk_reduction_lock
      if (rk_reduction_lock > 0) rk_reduction_lock = rk_reduction_lock - 1
      
      ! reset to the rank of the approximation which we use outside of the integrator
      X%rk = rk - 1      
   end subroutine rank_adaptive_PS_DLRA_riccati_step_rdp

   subroutine G_forward_map_riccati_rdp(X, B, CT, Qc, Rinv, tau, info, ifpred, T0, Tt, U0, Ut)
      !! This subroutine computes the solution of the non-stiff non-linear part of the 
      !! differential equation numerically using first-order explicit Euler.
      !! The update of the full low-rank factorization requires three separate
      !! steps called K, S, L.
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_vector_rdp),             intent(in)    :: B(:)
      !! System input.
      class(abstract_vector_rdp),             intent(in)    :: CT(:)
      !! System output.
      real(wp),                               intent(in)    :: Qc(:,:)
      !! Measurement weights.
      real(wp),                               intent(in)    :: Rinv(:,:)
      !! Inverse of the actuation weights.
      real(wp),                               intent(in)    :: tau
      !! Time step.
      integer,                                intent(out)   :: info
      !! Information flag.
      logical,                optional,       intent(in)    :: ifpred
      !! For Strang splitting: Determine whether we are in the predictor or corrector step
      class(abstract_vector_rdp), optional,   intent(inout) :: T0(:)  ! will be reused as Gamma
      class(abstract_vector_rdp), optional,   intent(in)    :: Tt(:)
      class(abstract_vector_rdp), optional,   intent(in)    :: U0(:)
      class(abstract_vector_rdp), optional,   intent(in)    :: Ut(:)
      !! Intermediate values
      
      ! Internal variables
      integer                                               :: rk

      if (time_lightROM()) call lr_timer%start('G_forward_map_riccati_rdp')

      rk = size(X%U)
      if (.not. allocated(U1))  allocate(U1( 1:rk), source=X%U(1)); 
      if (.not. allocated(QU))  allocate(QU( 1:rk), source=X%U(1));
      if (.not. allocated(Swrk0)) allocate(Swrk0(1:rk,1:rk))
      call zero_basis(U1); call zero_basis(QU); Swrk0 = 0.0_wp

      if (present(ifpred)) then
         ! second order in time
         if (ifpred) then ! predictor step with precomputed T0
            call K_step_riccati(X, U1, QU, B, CT, Qc, Rinv,     tau, info, reverse=.false., NL=T0)
            call S_step_riccati(X, U1, QU, B, CT, Qc, Rinv,     tau, info, reverse=.false.)
            call L_step_riccati(X, U1,     B, CT, Qc, Rinv,     tau, info)
         else             ! corrector step with precomputed T0, Tt and U0, Ut
            ! forward steps based on T0, U0 (within X)
            call K_step_riccati(X, U1, QU, B, CT, Qc, Rinv, 0.5*tau, info, reverse=.false., NL=T0)
            call S_step_riccati(X, U1, QU, B, CT, Qc, Rinv, 0.5*tau, info, reverse=.false.)
            call L_step_riccati(X, U1,     B, CT, Qc, Rinv,     tau, info)
            ! Compute Gamma = 0.5*(T0 @ (U1.T @ U0) + Tt @ (U1.T @ Ut))
            call zero_basis(QU); Swrk0 = 0.0_wp                         ! we use QU as a scratch array
            block
               class(abstract_vector_rdp), allocatable :: Xwrk(:)
               Swrk0 = innerprod(X%U, U0)
               call linear_combination(Xwrk, T0, Swrk0); call copy(QU, Xwrk)
               Swrk0 = innerprod(X%U, Ut)
               call linear_combination(Xwrk, Tt, Swrk0); call copy(T0, Xwrk) ! overwrite T0 with Gamma
            end block
            call axpby_basis(0.5_wp, QU, 0.5_wp, T0)
            ! Update X to most recent value
            call copy(X%U, U1)
            ! reverse steps based on Gamma
            call S_step_riccati(X, U1, QU, B, CT, Qc, Rinv, 0.5*tau, info, reverse=.true., NL=T0)
            call K_step_riccati(X, U1, QU, B, CT, Qc, Rinv, 0.5*tau, info, reverse=.true., NL=T0)
         end if
      else
         ! first order in time
         call K_step_riccati(   X, U1, QU, B, CT, Qc, Rinv,     tau, info)
         call S_step_riccati(   X, U1, QU, B, CT, Qc, Rinv,     tau, info)
         call L_step_riccati(   X, U1,     B, CT, Qc, Rinv,     tau, info)
      end if
      
      ! Copy updated low-rank factors to output
      call copy(X%U, U1)

      if (time_lightROM()) call lr_timer%stop('G_forward_map_riccati_rdp')
   end subroutine G_forward_map_riccati_rdp

   subroutine K_step_riccati_rdp(X, U1, QU, B, CT, Qc, Rinv, tau, info, reverse, NL)
      class(abstract_sym_low_rank_state_rdp),intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_vector_rdp),            intent(out)    :: U1(:)
      !! Intermediate low-rank factor.
      class(abstract_vector_rdp),            intent(inout)  :: QU(:)
      !! Precomputed application of the inhomogeneity.
      class(abstract_vector_rdp),            intent(in)     :: B(:)
      !! System input.
      class(abstract_vector_rdp),            intent(in)     :: CT(:)
      !! System output.
      real(wp),                              intent(in)     :: Qc(:,:)
      !! Measurement weights.
      real(wp),                              intent(in)     :: Rinv(:,:)
      !! Inverse of the actuation weights.
      real(wp),                              intent(in)     :: tau
      !! Time step.
      integer,                               intent(out)    :: info
      !! Information flag.
      logical, optional,                     intent(in)     :: reverse
      !! For Strang splitting: Determine if we are in forward or reverse branch
      class(abstract_vector_rdp), optional,  intent(in)     :: NL(:)
      !! Precomputed non-linear term.

      ! Internal variables
      integer                                               :: rk
      logical                                               :: reverse_order

      if (time_lightROM()) call lr_timer%start('K_step_riccati_rdp')

      ! Optional arguments
      reverse_order = optval(reverse, .false.)

      info = 0

      rk = size(X%U)
      if (.not. allocated(Uwrk0)) allocate(Uwrk0(1:rk), source=X%U(1)); 
      if (.not. allocated(Swrk0)) allocate(Swrk0(1:rk,1:rk)); 
      call zero_basis(Uwrk0); Swrk0 = 0.0_wp

      ! Constant part --> QU
      if (.not. reverse_order) then
         ! compute QU and pass to S step
         call apply_outerprod_w(QU, X%U, CT, Qc)
      end if

      ! Non-linear part --> Uwrk0
      block
         class(abstract_vector_rdp), allocatable :: Xwrk(:)
         call linear_combination(Xwrk, X%U, X%S)                       ! K0 = U0 @ S0
         call copy(U1, Xwrk)
      end block
      if (.not.present(NL)) then
         call apply_premult_outerprod_w(Swrk0, X%U, U1, B, Rinv)  ! (U0.T) @ B @ R^(-1) @ B.T @ K0
         block
            class(abstract_vector_rdp), allocatable :: Xwrk(:)
            call linear_combination(Xwrk, U1, Swrk0)   ! K0 @ Swrk0
            call copy(Uwrk0, Xwrk)
         end block
      else  ! non-linear term precomputed
         call copy(Uwrk0, NL)
      end if
      
      ! Combine to form G( K @ U.T ) @ U --> Uwrk0
      call axpby_basis(1.0_wp, QU, -1.0_wp, Uwrk0)

      ! Construct intermediate solution U1
      call axpby_basis(tau, Uwrk0, 1.0_wp, U1)                   ! K0 + tau*Kdot

      ! Orthonormalize in-place
      call qr(U1, Swrk0, info)
      call check_info(info, 'qr', this_module, 'K_step_Riccati_rdp')
      X%S = Swrk0

      if (time_lightROM()) call lr_timer%stop('K_step_riccati_rdp')
   end subroutine K_step_riccati_rdp

   subroutine S_step_riccati_rdp(X, U1, QU, B, CT, Qc, Rinv, tau, info, reverse, NL)
      class(abstract_sym_low_rank_state_rdp),intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_vector_rdp),            intent(in)    :: U1(:)
      !! Intermediate low-rank factor.
      class(abstract_vector_rdp),            intent(inout) :: QU(:)
      !! Precomputed application of the inhomogeneity.
      class(abstract_vector_rdp),            intent(in)    :: B(:)
      !! System input.
      class(abstract_vector_rdp),            intent(in)    :: CT(:)
      !! System output.
      real(wp),                              intent(in)    :: Qc(:,:)
      !! Measurement weights.
      real(wp),                              intent(in)    :: Rinv(:,:)
      !! Inverse of the actuation weights.
      real(wp),                              intent(in)    :: tau
      !! Time step.
      integer,                               intent(out)   :: info
      !! Information flag.
      logical, optional,                     intent(in)    :: reverse
      !! For Strang splitting: Determine if we are in forward or reverse branch
      class(abstract_vector_rdp), optional,  intent(in)    :: NL(:)
      !! Precomputed non-linear term.

      ! Internal variables
      integer                                              :: rk
      logical                                              :: reverse_order

      if (time_lightROM()) call lr_timer%start('S_step_riccati_rdp')

      info = 0

      ! Optional arguments
      reverse_order = optval(reverse, .false.)

      rk = size(X%U)
      if (.not. allocated(Swrk0)) allocate(Swrk0(1:rk,1:rk))
      if (.not. allocated(Swrk1)) allocate(Swrk1(1:rk,1:rk))
      Swrk0 = 0.0_wp; Swrk1 = 0.0_wp

      ! Constant part --> Swrk0
      if (reverse_order) then
         ! Compute QU and pass to K step
         call apply_outerprod_w(QU, X%U, CT, Qc)
      endif
      Swrk0 = innerprod(U1, QU)

      ! Non-linear part --> Swrk1
      if (.not.present(NL)) then
         call apply_premult_outerprod_w(Swrk1, X%U, U1, B, Rinv) !       U0.T @ B @ R^(-1) @ B.T @ U1
         Swrk1 = matmul(X%S, matmul(Swrk1, X%S))                 ! Sh @ (U0.T @ B @ R^(-1) @ B.T @ U1) @ Sh
      else ! Non-linear term precomputed
         Swrk1 = innerprod(U1, NL)
      end if

      ! Combine to form -U1.T @ G( U1 @ Sh @ U0.T ) @ U0
      Swrk0 = Swrk1 - Swrk0

      ! Construct intermediate coefficient matrix
      X%S = X%S + tau*Swrk0

      if (time_lightROM()) call lr_timer%stop('S_step_riccati_rdp')
   end subroutine S_step_riccati_rdp

   subroutine L_step_riccati_rdp(X, U1, B, CT, Qc, Rinv, tau, info)
      class(abstract_sym_low_rank_state_rdp),intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_vector_rdp),            intent(in)    :: U1(:)
      !! Intermediate low-rank factor.
      class(abstract_vector_rdp),            intent(in)    :: B(:)
      !! System input.
      class(abstract_vector_rdp),            intent(in)    :: CT(:)
      !! System output.
      real(wp),                              intent(in)    :: Qc(:,:)
      !! Measurement weights.
      real(wp),                              intent(in)    :: Rinv(:,:)
      !! Inverse of the actuation weights.
      real(wp),                              intent(in)    :: tau
      !! Time step.
      integer,                               intent(out)   :: info
      !! Information flag.

      ! Internal variables
      integer                                              :: rk

      if (time_lightROM()) call lr_timer%start('L_step_riccati_rdp')

      info = 0

      rk = size(X%U)
      if (.not. allocated(Uwrk0)) allocate(Uwrk0(1:rk), source=X%U(1))
      if (.not. allocated(Uwrk1)) allocate(Uwrk1(1:rk), source=X%U(1))
      if (.not. allocated(Swrk0)) allocate(Swrk0(1:rk,1:rk)) 
      call zero_basis(Uwrk0); call zero_basis(Uwrk1); Swrk0 = 0.0_wp

      block
         class(abstract_vector_rdp), allocatable :: Xwrk(:)
         call linear_combination(Xwrk, X%U, transpose(X%S))  ! L0.T                             U0 @ S.T
         call copy(Uwrk1, Xwrk)
      end block
      ! Constant part --> Uwrk0
      call apply_outerprod_w(Uwrk0, U1, CT, Qc)

      ! Non-linear part --> U
      call apply_premult_outerprod_w(Swrk0, U1, Uwrk1, B, Rinv)  !    U1.T @ B @ R^(-1) @ B.T @ U0 @ S.T
      block
         class(abstract_vector_rdp), allocatable :: Xwrk(:)
         call linear_combination(Xwrk, Uwrk1, Swrk0)  ! (U0 @ S.T) @ (U1.T @ B @ R^(-1) @ B.T @ U0 @ S.T)
         call copy(X%U, Xwrk)
      end block

      ! Combine to form U1.T @ G( U1.T@L.T )
      call axpby_basis(-1.0_wp, X%U, 1.0_wp, Uwrk0)

      ! Construct solution L1.T
      call axpby_basis(tau, Uwrk0, 1.0_wp, Uwrk1)               ! L0.T + tau*Ldot.T

      ! Update coefficient matrix
      X%S = innerprod(Uwrk1, U1)

      if (time_lightROM()) call lr_timer%stop('L_step_riccati_rdp')
   end subroutine L_step_riccati_rdp

   subroutine set_initial_rank_riccati(X, A, B, CT, Qc, Rinv, tau, mode, exptA, trans, tol, rk_init, nsteps)
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_linop_rdp),              intent(inout) :: A
      !! Linear operator
      class(abstract_vector_rdp),             intent(in)    :: B(:)
      !! System input.
      class(abstract_vector_rdp),             intent(in)    :: CT(:)
      !! System output.
      real(wp),                               intent(in)    :: Qc(:,:)
      !! Measurement weights.
      real(wp),                               intent(in)    :: Rinv(:,:)
      !! Inverse of the actuation weights.
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
      integer,                      optional, intent(in)    :: rk_init
      !! Smallest tested rank
      integer,                      optional, intent(in)    :: nsteps
      integer                                               :: n
      !! Number of steps to run before checking the singular values

      ! internal
      integer                                               :: i, irk, info, rkmax
      class(abstract_vector_rdp),               allocatable :: Utmp(:)
      real(wp),                                 allocatable :: Stmp(:,:), svals(:)
      logical                                               :: found, accept_rank
      character(len=512)                                    :: msg, fmt

      ! optional arguments
      X%rk = optval(rk_init, 1)
      n = optval(nsteps, 5)
      rkmax = size(X%U)

      info = 0
      accept_rank = .false.

      ! save initial condition
      allocate(Utmp(rkmax), source=X%U)
      allocate(Stmp(rkmax,rkmax)); Stmp = X%S

      do while (.not. accept_rank .and. X%rk <= rkmax)
         write(msg,'(4X,A,I0)') 'Test r = ', X%rk
         call log_message(msg, this_module, 'set_initial_rank_riccati')
         ! run integrator
         do i = 1,n
            call projector_splitting_DLRA_riccati_step_rdp(X, A, B, CT, Qc, Rinv, tau, mode, info, exptA, trans)
         end do

         ! check if singular values are resolved
         svals = svdvals(X%S(:X%rk,:X%rk))
         found = .false.
         tol_chk: do irk = 1, X%rk
            if ( svals(irk) < tol ) then
               found = .true.
               exit tol_chk
            end if
         end do tol_chk
         if (.not. found) irk = irk - 1
         write(msg,'(4X,A,I2,A,E8.2)') 'rk = ', X%rk, ' s_r =', svals(X%rk)
         call log_debug(msg, this_module, 'set_initial_rank_riccati')
         if (found) then
            accept_rank = .true.
            X%rk = irk
            write(msg,'(4X,A,I2,A,E10.4)') 'Accpeted rank: r = ', X%rk-1, ',     s_{r+1} = ', svals(X%rk)
            call log_message(msg, this_module, 'set_initial_rank_riccati')
         else
            X%rk = 2*X%rk
         end if
         
         ! reset initial conditions
         call copy(X%U, Utmp)
         X%S = Stmp
      end do

      if (X%rk > rkmax) then
         write(msg, *) 'Maximum rank reached but singular values are not converged. Increase rkmax and restart.'
         call stop_error(msg, this_module, 'set_initial_rank_riccati')
      end if

      ! reset to the rank of the approximation which we use outside of the integrator & mark rank as initialized
      X%rk = X%rk - 1
      X%rank_is_initialised = .true.
   end subroutine set_initial_rank_riccati

   subroutine reset_riccati_solver()
      ! internal
      character(len=128) :: msg
      write(msg,'(A,I0,A)') 'Riccati solver called ', RiccSolver_counter, ' times. Resetting coutner to 0.'
      call log_message(msg, this_module, 'DLRA_main')
      RiccSolver_counter = 0
      call reset_logfiles(logfile_basename)
   end subroutine reset_riccati_solver

end module LightROM_RiccatiSolvers
