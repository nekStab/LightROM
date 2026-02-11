module LightROM_DLRAIntegrators
   !! This module provides the implementation of the Krylov-based solvers for the Differential Lyapunov
   !! equation based on the dynamic low-rank approximation and operator splitting.
   ! Standard library
   use stdlib_linalg, only : eye, diag, svd, svdvals
   use stdlib_optval, only : optval
   ! LightKrylov modules
   use LightKrylov
   use LightKrylov, only: dp
   use LightKrylov_Logger
   use LightKrylov_AbstractVectors
   use LightKrylov_ExpmLib
   use LightKrylov_BaseKrylov
   ! LightROM modules
   use LightROM_AbstractLTIsystems
   use LightROM_SolverUtils
   use LightROM_LyapunovUtils
   use LightROM_Utils
   use LightROM_LoggerUtils
   use LightROM_Timing, only: lr_timer => global_lightROM_timer, time_lightROM
   
   implicit none
   character(len=*), parameter, private :: this_module = 'LR_LyapSolvers'

   character(len=*), parameter :: logfile_basename_lyapunov = 'Lyapunov_'
   integer :: LyapunovSolver_counter = 0
   character(len=*), parameter :: logfile_basename_riccati  = 'Riccati_'
   integer :: RiccatiSolver_counter = 0

   public :: Lyapunov_integrator
   public :: reset_solver

   interface Lyapunov_integrator
      module procedure Lyapunov_Integrator_rdp
   end interface

   interface Riccati_Integrator
      module procedure Riccati_Integrator_rdp
   end interface

contains

   subroutine Lyapunov_integrator_rdp(X, A, B, Tend, tau, info, exptA, iftrans, options)
      !! Main driver for the numerical integrator for the matrix-valued differential Lyapunov equation of the form
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
      !!   dynamical low-rank approximation to the Lyapunov equations by Mena et al. (2018). As the solutions
      !!   to the Lyapunov equation are by construction SPD, this fact can be taken advantage of to reduce the 
      !!   computational cost of the integration and, in particular, doing away with one QR factorization per timestep
      !!   while maintaining symmetry of the resulting matrix factorization.
      !! - The final element is the addition of the capability of dyanmic rank adaptivity for the projector-splitting
      !!   integrator proposed by Hochbruck et al. (2023). At the cost of integrating a supplementary solution vector, 
      !!   the rank of the solution is dynamically adapted to ensure that the corresponding additional singular value
      !!   stays below a chosen threshold.
      !!
      !! **Algorithmic Features**
      !! 
      !! - Separate integration of the stiff inhomogeneous part of the Lyapunov equation and the non-stiff inhomogeneity
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
      !! Linear operator
      class(abstract_vector_rdp),              intent(in)    :: B(:)
      !! Low-Rank inhomogeneity.
      real(dp),                                intent(in)    :: Tend
      !! Integration time horizon.
      real(dp),                                intent(in)    :: tau
      !! Desired time step. The avtual time-step will be computed such as to reach Tend in an integer number
      !! of steps.
      integer,                                 intent(out)   :: info
      !! Information flag
      procedure(abstract_exptA_rdp)                          :: exptA
      !! Routine for computation of the exponential propagator (default: Krylov-based exponential operator).
      logical,                       optional, intent(in)    :: iftrans
      !! Determine whether \(\mathbf{A}\) (default `.false.`) or \( \mathbf{A}^T\) (`.true.`) is used.
      type(dlra_opts),               optional, intent(in)    :: options
      !! Options for solver configuration
      
      ! internal
      character(len=*), parameter :: this_procedure = 'Lyapunov_integrator_rdp'
      character(len=*), parameter :: equation = 'Lyapunov'
      character(len=128) :: logfile
      type(dlra_opts) :: opts
      logical :: trans

      trans = optval(iftrans, .false.)
      logfile = logfile_basename_lyapunov
      LyapunovSolver_counter = LyapunovSolver_counter + 1

      call generic_initializer_rdp(X, opts, Tend, tau, rank_initializer, this_procedure, logfile, trans, options)

      call generic_DLRA_integrator_rdp(X, opts, Tend, tau, info, stepper, this_procedure, equation)

   contains

      subroutine stepper(if_rank_adaptive, info)
         logical, intent(in) :: if_rank_adaptive
         integer, intent(out) :: info
         !! Information flag
         if (if_rank_adaptive) then
            call rank_adaptive_projector_splitting_DLRA_step_lyapunov_rdp( &
                     &  X, A, B, opts%tau, opts%mode, exptA, trans, info, opts)
         else
            call fixed_rank_projector_splitting_DLRA_step_lyapunov_rdp( &
                     &  X, A, B, opts%tau, opts%mode, info, exptA, trans)
         end if
      end subroutine

      subroutine rank_initializer(X, opts)
         class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
         type(dlra_opts),                        intent(in)    :: opts
         ! internal
         integer :: rk_init, nsteps
         rk_init = 1
         nsteps = 5
         call set_initial_rank_lyapunov_rdp(X, A, B, opts%tau, opts%mode, exptA, trans, opts%tol, rk_init, nsteps)
      end subroutine

   end subroutine Lyapunov_integrator_rdp

   subroutine Riccati_Integrator_rdp(X, A, B, CT, Qc, Rinv, Tend, tau, info, exptA, iftrans, options)
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
      real(dp),                                intent(in)    :: Qc(:,:)
      !! Measurement weights.
      real(dp),                                intent(in)    :: Rinv(:,:)
      !! Inverse of the actuation weights.
      real(dp),                                intent(in)    :: Tend
      !! Integration time horizon. 
      real(dp),                                intent(in)    :: tau
      !! Desired time step. The avtual time-step will be computed such as to reach Tend in an integer number
      !! of steps.
      integer,                                 intent(out)   :: info
      !! Information flag.
      procedure(abstract_exptA_rdp)                          :: exptA
      !! Routine for computation of the exponential propagator (default: Krylov-based exponential operator).
      logical,                       optional, intent(in)    :: iftrans
      !! Determine whether \(\mathbf{A}\) (default `.false.`) or \( \mathbf{A}^T\) (`.true.`) is used.
      type(dlra_opts),               optional, intent(in)    :: options

      ! internal
      character(len=*), parameter :: this_procedure = 'Riccati_integrator_rdp'
      character(len=*), parameter :: equation = 'Riccati'
      character(len=128) :: logfile
      type(dlra_opts) :: opts
      logical :: trans

      trans = optval(iftrans, .false.)
      logfile = logfile_basename_riccati
      RiccatiSolver_counter = RiccatiSolver_counter + 1

      call generic_initializer_rdp(X, opts, Tend, tau, rank_initializer, this_procedure, logfile, trans, options)

      call generic_DLRA_integrator_rdp(X, opts, Tend, tau, info, stepper, this_procedure, equation)

   contains

      subroutine stepper(if_rank_adaptive, info)
         logical, intent(in) :: if_rank_adaptive
         integer, intent(out) :: info
         !! Information flag
         if (if_rank_adaptive) then
            call rank_adaptive_projector_splitting_DLRA_step_riccati_rdp( &
                     &  X, A, B, CT, Qc, Rinv, opts%tau, opts%mode, exptA, trans, opts, info)
         else
            call fixed_rank_projector_splitting_DLRA_step_riccati_rdp( &
                     &  X, A, B, CT, Qc, Rinv, opts%tau, opts%mode, info, exptA, trans)
         end if
      end subroutine

      subroutine rank_initializer(X, opts)
         class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
         type(dlra_opts),                        intent(in)    :: opts
         ! internal
         integer :: rk_init, nsteps
         rk_init = 1
         nsteps = 5
         call set_initial_rank_riccati_rdp(X, A, B, CT, Qc, Rinv, opts%tau, opts%mode, exptA, trans, opts%tol, rk_init, nsteps)
      end subroutine

   end subroutine Riccati_Integrator_rdp

   !!-------------------------------------------------------------------------------------------------------------
   !!
   !! Generic implementations of the initializer and integrator
   !!
   !!-------------------------------------------------------------------------------------------------------------

   subroutine generic_initializer_rdp(X, opts, Tend, tau, rank_initializer, this_procedure, logfile, trans, options)
      class(abstract_sym_low_rank_state_rdp),  intent(inout) :: X
      !! Low-Rank factors of the solution
      type(dlra_opts),                         intent(out)   :: opts
      !! Checked opts
      real(dp),                                intent(in)    :: Tend
      !! Integration time horizon.
      real(dp),                                intent(in)    :: tau
      !! Desired time step. The avtual time-step will be computed such as to reach Tend in an integer number
      !! of steps.
      procedure(abstract_dlra_rank_initializer)              :: rank_initializer
      !! Generic wrapper for the rank_initializer (Lyapnuov/Riccati)
      character(len=*),                        intent(in)    :: this_procedure
      !! Specific procedure name (for error attribution)
      character(len=*),                        intent(in)    :: logfile
      !! Basename for logfile
      logical,                                 intent(in)    :: trans
      !! Determine whether \(\mathbf{A}\) (default `.false.`) or \( \mathbf{A}^T\) (`.true.`) is used.
      type(dlra_opts),               optional, intent(in)    :: options
      !! Options for solver configuration

      ! internal
      character(len=128) :: msg
      integer :: rkmax

      call X%reset()
      ! Allocate memory for SVD & lagged fields
      rkmax = size(X%U)
      allocate(Usvd(rkmax,rkmax), ssvd(rkmax), VTsvd(rkmax,rkmax))

      ! Options
      if (present(options)) then
         opts = options
      else ! default
         opts = dlra_opts()
      end if
      opts%tau = tau
      opts%Tend = Tend
      call opts%init()

      call log_message('Initializing solver', this_module, this_procedure)
      write(msg,'(A,I0,A,F10.8)') 'Integration over ', opts%nsteps, ' steps with dt= ', opts%tau
      call log_information(msg, this_module, this_procedure)
      ! Prepare logfile
      call write_logfile_headers(logfile)

      ! determine initial rank if rank-adaptive
      if (opts%if_rank_adaptive .and. .not. X%rank_is_initialised) then
         call log_message('Determine initial rank:', this_module, this_procedure)
         call rank_initializer(X, opts)
      end if

      call log_settings(X, Tend, opts)
   end subroutine generic_initializer_rdp

   subroutine generic_DLRA_integrator_rdp(X, opts, Tend, tau, info,stepper, this_procedure, equation)
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      type(dlra_opts),                        intent(inout) :: opts
      !! DLRA opts
      real(dp),                               intent(in)    :: Tend
      !! Integration time horizon.
      real(dp),                               intent(in)    :: tau
      !! Desired time step. The avtual time-step will be computed such as to reach Tend in an integer number
      !! of steps.
      integer,                                intent(out)   :: info
      !! Information flag
      procedure(abstract_dlra_stepper)                      :: stepper
      !! Generic wrapper for the stepper (Lyapnuov/Riccati)
      character(len=*),                       intent(in)    :: this_procedure
      !! Specific procedure name (for error attribution)
      character(len=*),                       intent(in)    :: equation
      !! Solver type (Lyapunov/Riccati) - is used as in logfile name

      ! internal
      character(len=128) :: msg, logfile
      integer :: istep, irk, counter
      logical :: if_lastep
      real(dp), dimension(:), allocatable :: svals, svals_lag

      if (equation == 'Lyapunov') then
         logfile = logfile_basename_lyapunov
         counter = LyapunovSolver_counter
      else if (equation == 'Riccati') then
         logfile = logfile_basename_riccati
         counter = RiccatiSolver_counter
      else
         call stop_error('Unknown equation type: '//trim(equation), this_module, this_procedure)
      end if

      if (time_lightROM()) call lr_timer%start(this_procedure)
      call log_message('Starting DLRA integration', this_module, this_procedure)

      if_lastep = .false.
      if (opts%if_rank_adaptive) opts%rk_reduction_lock = opts%rk_reduction_barrier

      dlra : do istep = 1, opts%nsteps

         call log_step(X, istep, opts%nsteps)

         ! save lag data defore the timestep
         if ( X%rk > 0 .and. (mod(istep, opts%chkstep) == 0 .or. istep == opts%nsteps) ) then
            svals_lag = svdvals(X%S(:X%rk,:X%rk))
         end if

         ! dynamical low-rank approximation solver
         call stepper(opts%if_rank_adaptive, info)

         ! update time & step counters
         call X%increment_counters(opts%tau)

         ! here we can do some checks such as whether we have reached steady state
         if ( X%rk > 0 .and. (mod(istep, opts%chkstep) == 0 .or. istep == opts%nsteps) ) then
            svals = svdvals(X%S(:X%rk,:X%rk))
            if (.not. allocated(svals_lag)) allocate(svals_lag(X%rk), source=zero_rdp)
            call log_svals(logfile, X, opts%tau, svals, svals_lag, counter, istep, opts%nsteps)
            ! Check convergence
            if (istep == opts%nsteps) if_lastep = .true.
            irk = min(size(svals), size(svals_lag))
            X%is_converged = is_converged(X, svals(:irk), svals_lag(:irk), opts, if_lastep)
            if (X%is_converged) then
               write(msg,'(A,I0,A)') "Step ", istep, ": Solution converged!"
               call log_information(msg, this_module, this_procedure)
               exit dlra
            else ! if final step
               if (if_lastep) then
                  write(msg,'(A,I0,A)') "Step ", istep, ": Solution not converged!"
                  call log_information(msg, this_module, this_procedure)
               end if
            end if
         endif
      enddo dlra
      call log_message('Exiting solver', this_module, this_procedure)

   end subroutine generic_DLRA_integrator_rdp

end module LightROM_DLRAIntegrators
