module LightROM_LyapunovSolvers
   !! This module provides the implementation of the Krylov-based solvers for the Differential Lyapunov
   !! equation based on the dynamic low-rank approximation and operator splitting.
   ! Standard library
   use stdlib_linalg, only: eye, diag, svd, svdvals
   use stdlib_optval, only: optval
   use stdlib_logger, only: logger => global_logger
   ! LightKrylov modules
   use LightKrylov
   use LightKrylov, only: wp => dp
   use LightKrylov_Constants
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

   ! lagged solution for computation of increment norm
   class(abstract_vector_rdp),  allocatable   :: U_lag(:)
   real(wp),                    allocatable   :: S_lag(:,:)

   ! svd
   real(wp),                    allocatable   :: ssvd(:)

   ! module name
   private :: this_module
   character*128, parameter :: this_module = 'LightROM_LyapunovSolvers'

   public :: projector_splitting_DLRA_lyapunov_integrator
   public :: M_forward_map
   public :: G_forward_map_lyapunov
   public :: K_step_lyapunov
   public :: S_step_lyapunov
   public :: L_step_lyapunov

   interface projector_splitting_DLRA_lyapunov_integrator
      module procedure projector_splitting_DLRA_lyapunov_integrator_rdp
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

   subroutine projector_splitting_DLRA_lyapunov_integrator_rdp(X, A, B, Tend, tau, info, &
                                                                    & exptA, iftrans, options)
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
      real(wp),                                intent(in)    :: Tend
      !! Integration time horizon.
      real(wp),                                intent(inout) :: tau
      !! Desired time step. The avtual time-step will be computed such as to reach Tend in an integer number
      !! of steps.
      integer,                                 intent(out)   :: info
      !! Information flag
      procedure(abstract_exptA_rdp), optional                :: exptA
      !! Routine for computation of the exponential propagator (default: Krylov-based exponential operator).
      logical,                       optional, intent(in)    :: iftrans
      logical                                                :: trans
      !! Determine whether \(\mathbf{A}\) (default `.false.`) or \( \mathbf{A}^T\) (`.true.`) is used.
      type(dlra_opts),               optional, intent(in)    :: options
      type(dlra_opts)                                        :: opts
      !! Options for solver configuration

      ! Internal variables
      integer                                                :: istep, nsteps, rk
      integer                                                :: rk_reduction_lock   ! 'timer' to disable rank reduction
      real(wp)                                               :: tau0                ! input timestep
      real(wp)                                               :: nrm, nrmX, res      ! increment and solution norm
      real(wp)                                               :: El                  ! aggregate error estimate
      real(wp)                                               :: err_est             ! current error estimate
      real(wp)                                               :: tol                 ! current tolerance
      real(wp)                                               :: scale
      logical                                                :: verbose, converged
      character(len=4096)                                    :: msg
      character(len=128)                                     :: fmt
      procedure(abstract_exptA_rdp), pointer                 :: p_exptA => null()

      integer, allocatable :: log_units(:)
      integer :: i

      ! Optional arguments
      trans = optval(iftrans, .false.)

      ! Options
      if (present(options)) then
         opts = options
      else ! default
         opts = dlra_opts()
      end if

      if (present(exptA)) then
         p_exptA => exptA
      else
         p_exptA => k_exptA_rdp
      end if

      ! Compute number of steps
      tau0 = tau
      nsteps = nint(Tend/tau)
      tau = Tend/nsteps
      if (tau0 - tau > rtol_dp .and. opts%verbose) then
         write(msg, *) 'Reset timestep dt:', tau0, '-->', tau 
         if (io_rank()) call logger%log_message(trim(msg), module=this_module, procedure='RA-DLRA')
      end if

      ! Initialize
      rk_reduction_lock = 10
      converged         = .false.
      scale             = X%U(1)%get_size()
           
      ! check opts
      call read_opts(opts, tau, X%rank_is_initialized)

      ! set tolerance (might be overridden by error estimate)
      tol = opts%tol
      
      ! determine initial rank if rank-adaptive
      if (opts%if_rank_adaptive) then
         if (.not. X%rank_is_initialized) then
            call set_initial_rank(X, A, B, tau, opts%mode, p_exptA, trans, opts)
         else
            if (opts%verbose) then
               write(msg, '(A,F6.2,A,I4)') 'T = ', X%time,': Initial rank set to rk = ', X%rk
               if (io_rank()) call logger%log_message(trim(msg), module=this_module, procedure='RA-DLRA')
            end if
         end if
         if (opts%use_err_est) then
            El      = 0.0_wp
            err_est = splitting_error(X, A, B, tau, opts%mode, exptA, trans)
            tol = err_est / sqrt(X%U(1)%get_size() - real(X%rk + 1))
            if (opts%verbose) then
               write(msg, *) 'Initialization complete: rk = ', X%rk, ', local error estimate: ', tol
               if (io_rank()) call logger%log_message(trim(msg), module=this_module, procedure='RA-DLRA')
            end if
         end if
      else
         if (opts%verbose) then
            write(msg, *) 'Constant rank. rk =', X%rk
            if (io_rank()) call logger%log_message(trim(msg), module=this_module, procedure='RA-DLRA')
         end if
      end if

      if (opts%verbose) then
         write(msg, '(A)') 'DLRA initialization complete.'
         if (io_rank()) call logger%log_message(trim(msg), module=this_module, procedure='DLRA')
         write(msg, '(A,I6,A,E8.2)') 'Integration over ', nsteps, ' steps with dt = ', tau
         if (io_rank()) call logger%log_message(trim(msg), module=this_module, procedure='DLRA')
      end if

      dlra : do istep = 1, nsteps
         if ( opts%chkstep == 1 .or. nsteps == 1 ) then
            if (allocated(U_lag)) deallocate(U_lag)
            if (allocated(S_lag)) deallocate(S_lag)
            ! allocate lag data (we do it here so we do not need to store the data size and can pass the whole array)
            allocate(U_lag(X%rk), source=X%U(:X%rk)) ! U_lag = X%U
            allocate(S_lag(X%rk, X%rk)); S_lag = X%S(:X%rk,:X%rk)
            if (opts%verbose) then
               write(msg, *) 'Solution saved for increment norm computation.'
               if (io_rank()) call logger%log_debug(trim(msg), module=this_module, procedure='DLRA')
            end if
         end if
         ! dynamical low-rank approximation solver
         if (opts%if_rank_adaptive) then
            call rank_adaptive_PS_DLRA_lyapunov_step_rdp(X, A, B, tau, opts%mode, info, rk_reduction_lock, & 
                                                      & p_exptA, trans, verbose, tol)
            if ( opts%use_err_est ) then
               if ( mod(istep, opts%err_est_step) == 0 ) then
                  El = El + splitting_error(X, A, B, tau, opts%mode, exptA, trans)
                  tol = El / sqrt(X%U(1)%get_size() - real(X%rk + 1))
                  if (opts%verbose) then
                     write(msg, '(3X,I3,A,E8.2)') istep, ': recomputed error estimate: ', tol
                     if (io_rank()) call logger%log_message(trim(msg), module=this_module, procedure='RA-DLRA')
                  end if
               else
                  El = El + err_est
               end if
            end if
            !
         else
            call projector_splitting_DLRA_lyapunov_step_rdp(X, A, B, tau, opts%mode, info,  & 
                                                           & p_exptA, trans, opts%verbose)
         end if

         ! update time
         X%time = X%time + tau

         ! save lag data
         if ( opts%chkstep > 1 .and. (mod(istep + 1, opts%chkstep) == 0 .or. istep == nsteps-1) ) then
            if (allocated(U_lag)) deallocate(U_lag)
            if (allocated(S_lag)) deallocate(S_lag)
            ! allocate lag data (we do it here so we do not need to store the data size and can pass the whole array)
            allocate(U_lag(X%rk), source=X%U(:X%rk)) ! U_lag = X%U
            allocate(S_lag(X%rk, X%rk)); S_lag = X%S(:X%rk,:X%rk)
            if (opts%verbose) then
               write(msg, *) 'Solution saved for increment norm computation.'
               if (io_rank()) call logger%log_debug(trim(msg), module=this_module, procedure='DLRA')
            endif
         end if

         ! here we can do some checks such as whether we have reached steady state
         if ( mod(istep, opts%chkstep) == 0 .or. istep == nsteps ) then
            nrmX = dense_frobenius_norm(X%S(:X%rk,:X%rk), scale)
            nrm  = increment_norm(X%U(:X%rk), X%S(:X%rk,:X%rk), U_lag, S_lag, scale)
            res  = 0.0_wp !CALE_res_norm(X, A, B)
            if (opts%print_svals) then
               if (opts%if_rank_adaptive) then
                  rk = X%rk + 1
               else
                  rk = X%rk
               end if               
               if (.not.allocated( ssvd)) allocate( ssvd(rk)); ssvd = 0.0_wp
               ssvd(:rk) = svdvals(X%S(:rk,:rk))
               write(fmt,'(A,I4,A)') '(A,I6,A,F6.2,A,I3,A,',rk,'(1X,E10.4))'
               write(msg, fmt) "  Step ", istep, ", T = ", X%time, ', rk = ', X%rk, ': ', ssvd(:rk)
               if (io_rank()) call logger%log_message(trim(msg), module=this_module, procedure='RA-DLRA')
            else
               if (opts%if_rank_adaptive) then
                  write(msg, '(A,I6,A,F6.2,A,I4)') "  Step ", istep, ", T = ", X%time, ": rk = ", X%rk
                  if (io_rank()) call logger%log_message(trim(msg), module=this_module, procedure='RA-DLRA')
               end if
            end if
            write(msg, '(A,I6,A,F6.2,A,E10.4,A,E10.4,A,E10.4,A,E10.4)') "Step ", istep, ", T = ", X%time, &
                    ": dX/dt = ", nrm/tau, ' X = ', nrmX, ' (dX/dt)/X = ', nrm/tau/nrmX, ' res = ', res
            if (io_rank()) call logger%log_message(trim(msg), module=this_module, procedure='DLRA_INFO')
            ! Check convergence
            if (opts%chk_convergence) then
               converged = is_converged(nrm, nrmX, opts)
               if (converged) then
                  write(msg, '(A,I5,A)') "Step ", istep, ": Solution converged!"
                  if (io_rank()) call logger%log_message(trim(msg), module=this_module, procedure='DLRA')
                  exit dlra
               else ! if final step
                  if (istep == nsteps) then
                     write(msg, '(A,I5,A)') "Step ", istep, ": Solution not converged to tolerance!"
                     if (io_rank()) call logger%log_message(trim(msg), module=this_module, procedure='DLRA')
                  end if
               end if
            end if
         end if
      enddo dlra

      if (allocated(U1)) deallocate(U1)
      if (allocated(Uwrk)) deallocate(Uwrk)
      if (allocated(BBTU)) deallocate(BBTU)
      if (allocated(Swrk)) deallocate(Swrk)
      if (allocated(ssvd)) deallocate(ssvd)

      return
   end subroutine projector_splitting_DLRA_lyapunov_integrator_rdp

   !-----------------------
   !-----     PSI     -----
   !-----------------------

   subroutine projector_splitting_DLRA_lyapunov_step_rdp(X, A, B, tau, mode, info, exptA, trans, verbose)
      !! Driver for the time-stepper defining the splitting logic for each step of the the 
      !! projector-splitting integrator
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_linop_rdp),              intent(inout) :: A
      !! Linear operator
      class(abstract_vector_rdp),             intent(in)    :: B(:)
      !! Low-Rank inhomogeneity.
      real(wp),                               intent(in)    :: tau
      !! Time step.
      integer,                                intent(in)    :: mode
      !! TIme integration mode. Only 1st (Lie splitting - mode 1) and 2nd (Strang splitting - mode 2) 
      !! orders are implemented.
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
   end subroutine projector_splitting_DLRA_lyapunov_step_rdp

   !-----------------------------
   !
   !     RANK-ADAPTIVE PSI 
   !
   !-----------------------------

   subroutine rank_adaptive_PS_DLRA_lyapunov_step_rdp(X, A, B, tau, mode, info, rk_reduction_lock, exptA, trans, verbose, tol)
      !! Wrapper for projector_splitting_DLRA_lyapunov_step_rdp adding the logic for rank-adaptivity
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_linop_rdp),              intent(inout) :: A
      !! Linear operator
      class(abstract_vector_rdp),             intent(in)    :: B(:)
      !! Low-Rank inhomogeneity.
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
      logical,                                intent(in)    :: verbose
      !! Toggle verbosity
      real(wp),                               intent(in)    :: tol
      !! Tolerance for rank-adaptivity

      ! Internal variables
      integer                                               :: istep, rk, rk_new, irk, nsteps_
      logical                                               :: accept_step, found
      real(wp),                               allocatable   :: coef(:)
      real(wp)                                              :: norm
      character*128                                         :: msg

      integer, parameter :: max_step = 10
      integer, parameter :: rkmin = 2

      ! Allocate memory for SVD
      if (.not.allocated( ssvd)) allocate( ssvd(size(X%U)))

      ! ensure that we are integrating one more rank than we use for approximation
      X%rk = X%rk + 1
      rk = X%rk ! this is only to make the code more readable
      
      accept_step = .false.
      istep = 1
      do while ( .not. accept_step .and. istep < max_step )
         ! run a regular step
         call projector_splitting_DLRA_lyapunov_step_rdp(X, A, B, tau, mode, info, exptA, trans, verbose)
         ! compute singular values of X%S
         ssvd = 0.0_wp
         ssvd(:rk) = svdvals(X%S(:rk,:rk))
         found = .false.
         tol_chk: do irk = 1, rk
            if ( ssvd(irk) < tol ) then
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
               call stop_error(trim(msg), module=this_module, procedure='rank_adaptive_PS_DLRA_lyapunov_step_rdp')
            else
               write(msg,'(A, F8.3, A, I4)') 'T = ', X%time, ': increase to rk =', rk
               if (io_rank()) call logger%log_debug(trim(msg), module=this_module, procedure='RA-DLRA')
               
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
                                 & procedure='rank_adaptive_PS_DLRA_lyapunov_step_rdp')
               call X%U(rk)%scal(1.0_wp / X%U(rk)%norm())

               rk_reduction_lock = 10 ! avoid rank oscillations

            end if
         else ! the rank of the solution is sufficient
            accept_step = .true.

            if (irk /= rk .and. rk_reduction_lock == 0) then ! we should decrease the rank
               ! decrease rank

               rk_new = max(irk, rk - 2)  ! reduce by at most 2

               if (rk_new < rkmin) then
                  write(msg, '(A,I1,A,I1)') 'Cannot decrease rank, rkmin = ', rkmin, ' is reached. Resetting to previous rank: ', rk
                  if (io_rank()) call logger%log_warning(trim(msg), module=this_module, procedure='RA-DLRA')
               else
                  rk = rk_new
               end if

               write(msg, '(A, F8.3, A, I4)') 'T = ', X%time, ': decrease to rk =', rk - 1
               if (io_rank()) call logger%log_debug(trim(msg), module=this_module, procedure='RA-DLRA')
            end if
            
         end if ! found
         istep = istep + 1
      end do ! while .not. accept_step
      
      ! decrease rk_reduction_lock
      if (rk_reduction_lock > 0) rk_reduction_lock = rk_reduction_lock - 1
      
      ! reset to the rank of the approximation which we use outside of the integrator
      X%rk = rk - 1

      if (verbose) then
         write(msg,'(A,I4,A,I4,A,E14.8,A,I2)') 'rk = ', X%rk, ':     s_', irk,' = ', &
                                                   & ssvd(irk), ', rank_lock: ', rk_reduction_lock
         if (io_rank()) call logger%log_debug(trim(msg), module=this_module, procedure='RA-DLRA')
      end if

      if (.not. accept_step) then
         write(msg, *) 'Maximum number of substeps (', max_step, ') reached but new updated rank could not be found.'
         call stop_error(trim(msg), module=this_module, procedure='rank_adaptive_PS_DLRA_lyapunov_step_rdp')
      end if

      return
   end subroutine rank_adaptive_PS_DLRA_lyapunov_step_rdp

   subroutine set_initial_rank(X, A, B, tau, mode, exptA, trans, opts)
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
      type(dlra_opts),                        intent(inout) :: opts

      ! internal
      integer                                               :: i, irk, info, rkmax, nsteps
      class(abstract_vector_rdp),               allocatable :: Utmp(:)
      real(wp),                                 allocatable :: Stmp(:,:)
      logical                                               :: found, accept_rank
      character*1024                                         :: msg, fmt

      ! optional arguments
      X%rk = 1
      rkmax = size(X%U)            

      ! Allocate memory for SVD
      if (.not.allocated( ssvd)) allocate( ssvd(rkmax))

      info = 0
      accept_rank = .false.

      ! save initial condition
      allocate(Utmp(rkmax), source=X%U)
      allocate(Stmp(rkmax,rkmax)); Stmp = X%S
      
      do while (.not. accept_rank .and. X%rk <= rkmax-1)
         ! run integrator
         do i = 1, opts%ninit
            call projector_splitting_DLRA_lyapunov_step_rdp(X, A, B, tau, mode, info, exptA, trans, opts%verbose)
            if (opts%verbose) then
               ssvd = 0.0_wp
               ssvd(:X%rk) = svdvals(X%S(:X%rk,:X%rk))
               write(fmt, '(A,I4,A)') "(A,I3,1X,I2,1X,",X%rk,"(E10.2))"
               write(msg, fmt) 'Step', i, X%rk, ssvd(:X%rk)
               if (io_rank()) call logger%log_debug(trim(msg), module=this_module, procedure='set_initial_rank')
            end if
         end do

         ! check if singular values are resolved
         ssvd = 0.0_wp
         ssvd(:X%rk) = svdvals(X%S(:X%rk,:X%rk))
         found = .false.
         tol_chk: do irk = 1, X%rk
            if ( ssvd(irk) < opts%tol ) then
               found = .true.
               exit tol_chk
            end if
         end do tol_chk
         if (.not. found) irk = irk - 1
         if (found) then
            accept_rank = .true.
            write(msg,'(A,I4,A,I4,A,E10.4,A,I4)') 'rk = ', X%rk, ',     s(', irk, ') = ', ssvd(irk), &
                                          & ' : accepted. Approximation rank =', irk-1
            if (io_rank()) call logger%log_message(trim(msg), module=this_module, procedure='set_initial_rank')
            X%rk = irk
         else
            if (opts%verbose) then
               write(msg,'(A,I4,A,I4,A,E10.4,A)') 'rk = ', X%rk, ',     s(', X%rk, ') = ', ssvd(X%rk), ' : rejected.'
               if (io_rank()) call logger%log_message(trim(msg), module=this_module, procedure='set_initial_rank')
            end if
            X%rk = 2*X%rk
         end if
         
         ! reset initial conditions
         call copy_basis(X%U, Utmp)
         X%S = Stmp
      end do

      if (X%rk > rkmax-1) then
         write(msg, *) 'Maximum rank reached but singular values are not converged. Increase rkmax and restart.'
         call stop_error(trim(msg), module=this_module, procedure='set_initial_rank')
      end if

      ! reset to the rank of the approximation which we use outside of the integrator & mark rank as initialized
      X%rk = X%rk - 1
      X%rank_is_initialized = .true.

   end subroutine set_initial_rank

   subroutine M_forward_map_rdp(X, A, tau, info, exptA, iftrans)
      !! This subroutine computes the solution of the stiff linear part of the 
      !! differential equation exactly using the matrix exponential.
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
      integer                                               :: i, rk

      ! Optional argument
      trans = optval(iftrans, .false.)

      rk = X%rk
      allocate(R(rk,rk)); R = 0.0_wp 
      allocate(perm(rk)); perm = 0

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
      X%S(:rk,:rk) = matmul(R, matmul(X%S(:rk,:rk), transpose(R)))

      return
   end subroutine M_forward_map_rdp

   subroutine G_forward_map_lyapunov_rdp(X, B, tau, info)
      !! This subroutine computes the solution of the non-stiff part of the 
      !! differential equation numerically using first-order explicit Euler.
      !! The update of the full low-rank factorization requires three separate
      !! steps called K, S, L.
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
      rkmax = size(X%U)
      if (.not. allocated(Uwrk)) allocate(Uwrk(rkmax), source=X%U(1))
      call zero_basis(Uwrk)

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

   real(wp) function splitting_error(X, A, B, tau, mode, exptA, trans, scale) result(err_est)
      !! This function estimates the splitting error of the integrator as a function of the chosen timestep.
      !! This error estimation can be integrated over time to give an estimate of the compound error due to 
      !! the splitting approach.
      !! This error can be used as a tolerance for the rank-adaptivity to ensure that the low-rank truncation 
      !! error is smaller than the splitting error.
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
      real(wp),                     optional, intent(in)    :: scale
      real(wp)                                              :: scale_
      !! Scaling factor
      
      ! internals
      ! save current state to reset it later
      class(abstract_vector_rdp),               allocatable :: Utmp(:)
      real(wp),                                 allocatable :: Stmp(:,:)
      ! first solution to compute the difference against
      class(abstract_vector_rdp),               allocatable :: U1(:)
      real(wp),                                 allocatable :: S1(:,:)
      ! projected bases
      real(wp),                                 allocatable :: V1(:,:), V2(:,:)
      ! projected difference
      real(wp),                                 allocatable :: D(:,:)
      integer                                               :: rx, r, info

      rx = X%rk
      r  = 2*rx

      ! save curret state
      allocate(Utmp(rx), source=X%U(:rx))
      allocate(Stmp(rx,rx)); Stmp = X%S(:rx,:rx)

      ! tau step
      call projector_splitting_DLRA_lyapunov_step_rdp(X, A, B, tau, mode, info, exptA, trans, verbose=.false.)
      ! save result
      allocate(U1(rx), source=X%U(:rx))
      allocate(S1(rx,rx)); S1 = X%S(:rx,:rx)

      ! reset curret state
      call copy_basis(X%U(:rx), Utmp)
      X%S(:rx,:rx) = Stmp

      ! tau/2 steps
      call projector_splitting_DLRA_lyapunov_step_rdp(X, A, B, 0.5*tau, mode, info, exptA, trans, verbose=.false.)
      call projector_splitting_DLRA_lyapunov_step_rdp(X, A, B, 0.5*tau, mode, info, exptA, trans, verbose=.false.)
      
      ! compute local error based on frobenius norm of difference
      err_est = 2**mode / (2**mode - 1) * increment_norm(X%U(:rx), X%S(:rx,:rx), U1(:rx), S1, scale_)

   end function splitting_error

end module LightROM_LyapunovSolvers
