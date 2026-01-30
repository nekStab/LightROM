module LightROM_LyapunovSolvers
   !! This module provides the implementation of the Krylov-based solvers for the Differential Lyapunov
   !! equation based on the dynamic low-rank approximation and operator splitting.
   ! Standard library
   use stdlib_linalg, only : eye, diag, svd, svdvals
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
   use LightROM_LoggerUtils
   use LightROM_Timing, only: lr_timer => global_lightROM_timer, time_lightROM
   
   implicit none

   ! global scratch arrays
   real(wp),                    allocatable   :: ssvd(:)
   real(wp),                    allocatable   :: Usvd(:,:), VTsvd(:,:)

   ! module name
   private :: this_module
   character(len=*), parameter :: this_module = 'LR_LyapSolvers'
   character(len=*), parameter :: logfile_basename = 'Lyap_SVD'
   integer :: LyapSolver_counter = 0

   public :: projector_splitting_DLRA_lyapunov_integrator
   public :: M_forward_map
   public :: G_forward_map_lyapunov
   public :: K_step_lyapunov
   public :: S_step_lyapunov
   public :: L_step_lyapunov
   public :: reset_lyapunov_solver

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

   subroutine projector_splitting_DLRA_lyapunov_integrator_rdp(X, A, B, Tend, tau, info, exptA, iftrans, options)
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

      if (time_lightROM()) call lr_timer%start('DLRA_lyapunov_integrator_rdp')
      LyapSolver_counter = LyapSolver_counter + 1

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
      X%is_converged = .false.
      X%time = 0.0_wp
      X%step = 0
      ! Compute number of steps
      if_lastep = .false.
      nsteps = nint(Tend/tau)

      rkmax = size(X%U)
      ! Allocate memory for SVD & lagged fields
      allocate(Usvd(rkmax,rkmax), ssvd(rkmax), VTsvd(rkmax,rkmax))

      call log_message('Initializing Lyapunov solver', this_module, 'DLRA_main')
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
      endif

      ! determine initial rank if rank-adaptive
      if (opts%if_rank_adaptive) then
         rk_reduction_lock = 10  ! initialize rank reduction lock
         if (.not. X%rank_is_initialised) then
            call log_message('Determine initial rank:', this_module, 'DLRA_main')
            call set_initial_rank_lyapunov(X, A, B, tau, opts%mode, exptA, trans, tol)
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
            call rank_adaptive_PS_DLRA_lyapunov_step_rdp(X, A, B, tau, opts%mode, info, rk_reduction_lock, & 
                                                         & exptA, trans, tol)
         else
            call projector_splitting_DLRA_lyapunov_step_rdp(X, A, B, tau, opts%mode, info, exptA, trans)
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
            call log_svals(logfile_basename, X, tau, svals, svals_lag, LyapSolver_counter, istep, nsteps)
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
      call log_message('Exiting Lyapunov solver', this_module, 'DLRA_main')
      ! Clean up scratch space
      deallocate(Usvd, ssvd, VTsvd)
      if (time_lightROM()) call lr_timer%stop('DLRA_lyapunov_integrator_rdp')
   end subroutine projector_splitting_DLRA_lyapunov_integrator_rdp

   !-----------------------
   !-----     PSI     -----
   !-----------------------

   subroutine projector_splitting_DLRA_lyapunov_step_rdp(X, A, B, tau, mode, info, exptA, trans)
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
      
      ! Internal variables
      integer                                               :: istep, nsteps
      character(len=128)                                    :: msg

      if (time_lightROM()) call lr_timer%start('DLRA_lyapunov_step_rdp')

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

      if (time_lightROM()) call lr_timer%stop('DLRA_lyapunov_step_rdp')
   end subroutine projector_splitting_DLRA_lyapunov_step_rdp

   !-----------------------------
   !
   !     RANK-ADAPTIVE PSI 
   !
   !-----------------------------

   subroutine rank_adaptive_PS_DLRA_lyapunov_step_rdp(X, A, B, tau, mode, info, rk_reduction_lock, exptA, trans, tol)
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
         call projector_splitting_DLRA_lyapunov_step_rdp(X, A, B, tau, mode, info, exptA, trans)
         ! compute singular values of X%S
         call svd(X%S(:rk,:rk), ssvd(:rk), Usvd(:rk,:rk), VTsvd(:rk,:rk))
         call find_rank(found, irk, ssvd, tol)
         
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
               call decrease_rank(X, Usvd, ssvd, rk)
               rk = max(irk, rk - 2)  ! reduce by at most 2
               write(msg, '(A,I0)') 'Rank decreased to rk= ', rk
               call log_message(msg, this_module, 'DLRA_main')
            end if
            
         end if ! found
         istep = istep + 1
      end do ! while .not. accept_step

      if (istep >= max_step) then
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
   end subroutine rank_adaptive_PS_DLRA_lyapunov_step_rdp

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
      integer                                               :: i, rk

      if (time_lightROM()) call lr_timer%start('M_forward_map_rdp')
      
      ! Optional argument
      trans = optval(iftrans, .false.)

      rk = X%rk
      allocate(R(rk,rk)); R = 0.0_wp

      ! Apply propagator to initial basis
      allocate(exptAU, source=X%U(1)); call exptAU%zero()
      do i = 1, rk
         call exptA(exptAU, A, X%U(i), tau, info, trans)
         call copy(X%U(i), exptAU) ! overwrite old solution
      end do
      ! Reorthonormalize in-place
      call qr(X%U(:rk), R, info)
      call check_info(info, 'qr', this_module, 'M_forward_map_rdp')
   
      ! Update coefficient matrix
      X%S(:rk,:rk) = matmul(R, matmul(X%S(:rk,:rk), transpose(R)))

      if (time_lightROM()) call lr_timer%stop('M_forward_map_rdp')
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
      class(abstract_vector_rdp),  allocatable              :: U1(:)
      class(abstract_vector_rdp),  allocatable              :: BBTU(:)
      integer                                               :: rk

      if (time_lightROM()) call lr_timer%start('G_forward_map_lyapunov_rdp')
    
      rk = X%rk
      allocate(  U1(rk), source=X%U(1)); call zero_basis(U1)
      allocate(BBTU(rk), source=X%U(1)); call zero_basis(BBTU)

      call K_step_lyapunov(X, U1, BBTU, B, tau, info)
      call S_step_lyapunov(X, U1, BBTU,    tau, info)
      call L_step_lyapunov(X, U1,       B, tau, info)
      
      ! Copy updated low-rank factors to output
      call copy(X%U(:rk), U1)

      if (time_lightROM()) call lr_timer%stop('G_forward_map_lyapunov_rdp')
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
      class(abstract_vector_rdp), allocatable :: Uwrk(:)
      integer                                 :: rk

      if (time_lightROM()) call lr_timer%start('K_step_lyapunov_rdp')

      rk = X%rk
      call linear_combination(Uwrk, X%U(:rk), X%S(:rk,:rk))  ! K0
      call copy(U1, Uwrk)
      call apply_outerprod(BBTU, B, X%U(:rk))                ! Kdot
      ! Construct intermediate solution U1
      call axpby_basis(tau, BBTU, 1.0_wp, U1)                ! K0 + tau*Kdot
      ! Orthonormalize in-place
      call qr(U1, X%S(:rk,:rk), info)
      call check_info(info, 'qr', this_module, 'K_step_lyapunov_rdp')

      if (time_lightROM()) call lr_timer%stop('K_step_lyapunov_rdp')
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
      integer                                               :: rk
      real(wp),                               allocatable   :: Swrk(:,:)

      if (time_lightROM()) call lr_timer%start('S_step_lyapunov_rdp')

      rk = X%rk
      allocate(Swrk(rk,rk)); Swrk = 0.0_wp
      Swrk = innerprod(U1, BBTU)          ! - Sdot
      ! Construct intermediate coefficient matrix
      X%S(:rk,:rk) = X%S(:rk,:rk) - tau*Swrk

      if (time_lightROM()) call lr_timer%stop('S_step_lyapunov_rdp')
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
      integer                                               :: rk
      class(abstract_vector_rdp),             allocatable   :: Uwrk(:)

      if (time_lightROM()) call lr_timer%start('L_step_lyapunov_rdp')

      rk = X%rk
      call linear_combination(Uwrk, X%U(:rk), transpose(X%S(:rk,:rk)))  ! L0.T
      ! Construct derivative
      call apply_outerprod(X%U(:rk), B, U1)       ! Ldot.T
      ! Construct solution L1.T
      call axpby_basis(tau, X%U(:rk), 1.0_wp, Uwrk)
      ! Update coefficient matrix
      X%S(:rk,:rk) = innerprod(Uwrk, U1)

      if (time_lightROM()) call lr_timer%stop('L_step_lyapunov_rdp')
   end subroutine L_step_lyapunov_rdp

   subroutine set_initial_rank_lyapunov(X, A, B, tau, mode, exptA, trans, tol, rk_init, nsteps)
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
      integer,                      optional, intent(in)    :: rk_init
      !! Smallest tested rank
      integer,                      optional, intent(in)    :: nsteps
      integer                                               :: n
      !! Number of steps to run before checking the singular values

      ! internal
      integer                                               :: i, rk, irk, info, rkmax
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
         call log_message(msg, this_module, 'set_initial_rank_lyapunov')
         ! run integrator
         do i = 1,n
            call projector_splitting_DLRA_lyapunov_step_rdp(X, A, B, tau, mode, info, exptA, trans)
         end do

         ! check if singular values are resolved
         svals = svdvals(X%S(:X%rk,:X%rk))
         found = .false.
         tol_chk: do irk = 1, rk
            if ( ssvd(irk) < tol ) then
               found = .true.
               exit tol_chk
            end if
         end do tol_chk
         if (.not. found) irk = irk - 1

         write(msg,'(4X,A,I2,A,E8.2)') 'rk = ', X%rk, ' s_r =', svals(X%rk)
         call log_debug(msg, this_module, 'set_initial_rank_lyapunov')
         if (found) then
            accept_rank = .true.
            X%rk = irk
            write(msg,'(4X,A,I2,A,E10.4)') 'Accpeted rank: r = ', X%rk-1, ',     s_{r+1} = ', svals(X%rk)
            call log_message(msg, this_module, 'set_initial_rank_lyapunov')
         else
            X%rk = 2*X%rk
         end if
         
         ! reset initial conditions
         call copy(X%U, Utmp)
         X%S = Stmp
      end do

      if (X%rk > rkmax) then
         write(msg, *) 'Maximum rank reached but singular values are not converged. Increase rkmax and restart.'
         call stop_error(msg, this_module, 'set_initial_rank_lyapunov')
      end if

      ! reset to the rank of the approximation which we use outside of the integrator & mark rank as initialized
      X%rk = X%rk - 1
      X%rank_is_initialised = .true.
   end subroutine set_initial_rank_lyapunov

   subroutine compute_splitting_error(err_est, X, A, B, tau, mode, exptA, trans)
      !! This function estimates the splitting error of the integrator as a function of the chosen timestep.
      !! This error estimation can be integrated over time to give an estimate of the compound error due to 
      !! the splitting approach.
      !! This error can be used as a tolerance for the rank-adaptivity to ensure that the low-rank truncation 
      !! error is smaller than the splitting error.
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
      call projector_splitting_DLRA_lyapunov_step_rdp(X, A, B, tau, mode, info, exptA, trans)
      ! save result
      allocate(U1(rx), source=X%U(:rx))
      allocate(S1(rx,rx)); S1 = X%S(:rx,:rx)

      ! reset curret state
      call copy(X%U(:rx), Utmp)
      X%S(:rx,:rx) = Stmp

      ! tau/2 steps
      call projector_splitting_DLRA_lyapunov_step_rdp(X, A, B, 0.5*tau, mode, info, exptA, trans)
      call projector_splitting_DLRA_lyapunov_step_rdp(X, A, B, 0.5*tau, mode, info, exptA, trans)

      ! compute common basis
      call project_onto_common_basis(V1, V2, U1(:rx), X%U(:rx))

      ! project second low-rank state onto common basis and construct difference
      allocate(D(r,r)); D = 0.0_wp
      D(    :rx,     :rx) = S1 - matmul(V1, matmul(X%S(:rx,:rx), transpose(V1)))
      D(rx+1:r ,     :rx) =    - matmul(V2, matmul(X%S(:rx,:rx), transpose(V1)))
      D(    :rx, rx+1:r ) =    - matmul(V1, matmul(X%S(:rx,:rx), transpose(V2)))
      D(rx+1:r , rx+1:r ) =    - matmul(V2, matmul(X%S(:rx,:rx), transpose(V2)))

      ! compute local error based on frobenius norm of difference
      err_est = 2**mode / (2**mode - 1) * sqrt( sum( svdvals(D) ** 2 ) )

      ! reset curret state
      call copy(X%U(:rx), Utmp)
      X%S(:rx,:rx) = Stmp
   end subroutine compute_splitting_error

   subroutine reset_lyapunov_solver()
      ! internal
      character(len=128) :: msg
      write(msg,'(A,I0,A)') 'Lyapunov solver called ', LyapSolver_counter, ' times. Resetting coutner to 0.'
      call log_message(msg, this_module, 'DLRA_main')
      LyapSolver_counter = 0
      call reset_logfiles(logfile_basename)
   end subroutine reset_lyapunov_solver

end module LightROM_LyapunovSolvers
