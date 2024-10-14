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
   
   implicit none

   ! global scratch arrays
   !class(abstract_vector_rdp),  allocatable   :: U1(:)
   !class(abstract_vector_rdp),  allocatable   :: Uwrk(:)
   !class(abstract_vector_rdp),  allocatable   :: BBTU(:)
   !real(wp),                    allocatable   :: Swrk(:,:)

   ! lagged solution for computation of increment norm
   class(abstract_vector_rdp),  allocatable   :: U_lag(:)
   real(wp),                    allocatable   :: S_lag(:,:)

   ! svd
   !real(wp),                    allocatable   :: ssvd(:)
   !real(wp),                    allocatable   :: Usvd(:,:), VTsvd(:,:)

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
      procedure(abstract_exptA_rdp), optional                :: exptA
      !! Routine for computation of the exponential propagator (default: Krylov-based exponential operator).
      logical,                       optional, intent(in)    :: iftrans
      logical                                                :: trans
      !! Determine whether \(\mathbf{A}\) (default `.false.`) or \( \mathbf{A}^T\) (`.true.`) is used.
      type(dlra_opts),               optional, intent(in)    :: options
      type(dlra_opts)                                        :: opts
      !! Options for solver configuration

      ! Internal variables
      integer                                                :: istep, nsteps, chkstep
      integer                                                :: rk_reduction_lock   ! 'timer' to disable rank reduction
      real(wp)                                               :: T                   ! simulation time
      real(wp)                                               :: nrm, nrmX           ! increment and solution norm
      real(wp)                                               :: El                  ! aggregate error estimate
      real(wp)                                               :: err_est             ! current error estimate
      real(wp)                                               :: tol                 ! current tolerance
      logical                                                :: converged
      character*128                                          :: msg, fmt
      procedure(abstract_exptA_rdp), pointer                 :: p_exptA => null()

      ! Optional arguments
      trans = optval(iftrans, .false.)

      ! Options
      if (present(options)) then
         opts = options
      else ! default
         opts = dlra_opts()
      end if

      ! set tolerance
      tol     = opts%tol

      if (present(exptA)) then
         p_exptA => exptA
      else
         p_exptA => k_exptA_rdp
      end if

      ! Initialize
      T                 = 0.0_wp
      rk_reduction_lock = 10
      converged         = .false.

      ! Compute number of steps
      nsteps = nint(Tend/tau)
      write(msg,'(A,I0,A,F10.8)') 'Integration over ', nsteps, ' steps with dt= ', tau
      call logger%log_information(msg, module=this_module, procedure='DLRA')
      ! Pretty output
      write(fmt,'(A,I0,A)') '(A,I', max(5,int(log10(real(nsteps)))+1) , ',A,F10.4,3(A,E12.5))'
     
      ! Determine IO step
      chkstep = get_chkstep(opts, tau)

      if ( opts%mode > 2 ) then
         write(msg,'(A)') "Time-integration order for the operator splitting of d > 2 &
                      & requires adjoint solves and is not implemented. Resetting torder = 2." 
         call logger%log_message(msg, module=this_module, procedure='DLRA')
      else if ( opts%mode < 1 ) then
         write(msg,'(A,I0)') "Invalid time-integration order specified: ", opts%mode
         call stop_error(msg, module=this_module, procedure='DLRA')
      endif

      ! determine initial rank if rank-adaptive
      if (opts%if_rank_adaptive) then
         if (.not. X%rank_is_initialised) then
            call set_initial_rank(X, A, B, tau, opts%mode, p_exptA, trans, tol)
         end if
         if (opts%use_err_est) then
            err_est = 0.0_wp
            El      = 0.0_wp
            call compute_splitting_error(err_est, X, A, B, tau, opts%mode, exptA, trans)
            tol = err_est / sqrt(X%U(1)%get_size() - real(X%rk + 1))
            write(msg, *) 'Initialization complete: rk = ', X%rk, ', local error estimate: ', tol
            call logger%log_information(msg, module=this_module, procedure='RA-DLRA')
         end if
      end if

      dlra : do istep = 1, nsteps

         ! save lag data defore the timestep
         if (mod(istep, chkstep) == 0 .or. istep == nsteps ) then
            ! allocate lag data (we do it here so we do not need to store the data size and can pass the whole array)
            if (allocated(U_lag)) deallocate(U_lag)
            if (allocated(S_lag)) deallocate(S_lag)
            allocate(U_lag(X%rk), source=X%U(:X%rk)) ! U_lag = X%U
            allocate(S_lag(X%rk, X%rk)); S_lag = X%S(:X%rk, :X%rk)
            write(msg,'(A,I0,A)') 'Step ', istep, ': Solution saved before timestep for increment norm computation.'
            call logger%log_debug(msg, module=this_module, procedure='DLRA')
         end if

         ! dynamical low-rank approximation solver
         if (opts%if_rank_adaptive) then
            call rank_adaptive_PS_DLRA_lyapunov_step_rdp(X, A, B, tau, opts%mode, info, rk_reduction_lock, & 
                                                      & p_exptA, trans, tol)  
            if ( opts%use_err_est ) then
               if ( mod(istep, opts%err_est_step) == 0 ) then
                  call compute_splitting_error(err_est, X, A, B, tau, opts%mode, exptA, trans)
                  El = El + err_est
                  tol = El / sqrt(256_wp - real(X%rk + 1))
                  write(msg,'(3X,I3,A,E8.2)') istep, ': recomputed error estimate: ', tol
                  call logger%log_information(msg, module=this_module, procedure='RA-DLRA')
               else
                  El = El + err_est
               end if
            end if
         else
            call projector_splitting_DLRA_lyapunov_step_rdp(X, A, B, tau, opts%mode, info, p_exptA, trans)
         end if

         ! update time
         T = T + tau

         ! here we can do some checks such as whether we have reached steady state
         if (mod(istep, chkstep) == 0 .or. istep == nsteps) then
            call compute_increment_norm(nrm, X%U(:X%rk), X%S(:X%rk,:X%rk), U_lag, S_lag, .true.)
            nrmX = compute_norm(X, .true.)
            write(msg,fmt) "Step ", istep, ": T = ", T, ": dX= ", nrm, ' X= ', nrmX, ' dX/dt/X= ', nrm/tau/nrmX
            call logger%log_information(msg, module=this_module, procedure='DLRA')
            deallocate(U_lag); deallocate(S_lag)
            ! Check convergence
            converged = is_converged(nrm, nrmX, opts)
            if (converged) then
               write(msg,'(A,I0,A)') "Step ", istep, ": Solution converged!"
               call logger%log_information(msg, module=this_module, procedure='DLRA')
               exit dlra
            else ! if final step
               if (istep == nsteps) then
                  write(msg,'(A,I0,A)') "Step ", istep, ": Solution not converged!"
                  call logger%log_information(msg, module=this_module, procedure='DLRA')
               end if
            end if
         endif
      enddo dlra

      ! Clean up scratch space
      !if (allocated(Swrk)) deallocate(Swrk)
      !if (allocated(  U1)) deallocate(  U1)
      !if (allocated(Uwrk)) deallocate(Uwrk)
      !if (allocated(BBTU)) deallocate(BBTU)

      return
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
      integer                                               :: istep, rk, irk
      logical                                               :: accept_step, found
      real(wp),                               allocatable   :: coef(:)
      real(wp)                                              :: norm
      character*128                                         :: msg

      integer, parameter                                    :: max_step = 5  ! might not be needed

      real(wp),                    allocatable   :: ssvd(:)
      real(wp),                    allocatable   :: Usvd(:,:), VTsvd(:,:)

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
         call projector_splitting_DLRA_lyapunov_step_rdp(X, A, B, tau, mode, info, exptA, trans)
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
         if (.not. found) irk = irk - 1
         
         ! choose action
         if (.not. found) then ! none of the singular values is below tolerance
            ! increase rank and run another step
            if (rk == size(X%U)) then ! cannot increase rank without reallocating X%U and X%S
               write(msg,'(A)') 'Cannot increase rank, rkmax is reached. Increase rkmax and restart!'
               call stop_error(msg, module=this_module, procedure='rank_adaptive_PS_DLRA_lyapunov_step_rdp')
            else
               write(msg,'(A,I0)') 'rk= ', rk + 1
               call logger%log_warning(msg, module=this_module, procedure='RA-DLRA')
               
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
               
               ! rotate basis onto principal axes
               block
                  class(abstract_vector_rdp), allocatable :: Xwrk(:)
                  call linear_combination(Xwrk, X%U(:rk), Usvd(:rk,:rk))
                  call copy_basis(X%U(:rk), Xwrk)
               end block
               X%S(:rk,:rk) = diag(ssvd(:rk))

               rk = max(irk, rk - 2)  ! reduce by at most 2

               write(msg, '(A,I0)') 'rk= ', rk
               call logger%log_warning(msg, module=this_module, procedure='RA-DLRA')
            end if
            
         end if ! found
         istep = istep + 1
      end do ! while .not. accept_step
      write(msg,'(A,I3,A,I2,A,E14.8,A,I2)') 'rk = ', X%rk-1, ':     s_', irk,' = ', &
                                                   & ssvd(irk), ', rank_lock: ', rk_reduction_lock
      call logger%log_information(msg, module=this_module, procedure='RA-DLRA')

      ! decrease rk_reduction_lock
      if (rk_reduction_lock > 0) rk_reduction_lock = rk_reduction_lock - 1
      
      ! reset to the rank of the approximation which we use outside of the integrator
      X%rk = rk - 1      

      return
   end subroutine rank_adaptive_PS_DLRA_lyapunov_step_rdp

   subroutine set_initial_rank(X, A, B, tau, mode, exptA, trans, tol, rk_init, nsteps)
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
         svals = svdvals(X%S(:X%rk,:X%rk))
         ! run integrator
         do i = 1,n
            call projector_splitting_DLRA_lyapunov_step_rdp(X, A, B, tau, mode, info, exptA, trans)
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
         write(msg,'(A,I2,A,E8.2)') '   rk = ', X%rk, ' s_r =', svals(X%rk)
         call logger%log_debug(msg, module=this_module, procedure='set_initial_rank')
         if (found) then
            accept_rank = .true.
            X%rk = irk
            write(msg,'(A,I2,A,E10.4)') '   Accpeted rank: r = ', X%rk-1, ',     s_{r+1} = ', svals(X%rk)
            call logger%log_information(msg, module=this_module, procedure='set_initial_rank')
         else
            X%rk = 2*X%rk
         end if
         
         ! reset initial conditions
         call copy_basis(X%U, Utmp)
         X%S = Stmp
      end do

      if (X%rk > rkmax) then
         write(msg, *) 'Maximum rank reached but singular values are not converged. Increase rkmax and restart.'
         call stop_error(msg, module=this_module, procedure='set_initial_rank')
      end if

      ! reset to the rank of the approximation which we use outside of the integrator & mark rank as initialized
      X%rk = X%rk - 1
      X%rank_is_initialised = .true.
      
      return
   end subroutine set_initial_rank

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

      ! svd
      real(wp), dimension(:),                   allocatable :: ssvd

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
      call copy_basis(X%U(:rx), Utmp)
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
      
      ! svd
      ssvd = svdvals(D)

      ! compute local error based on frobenius norm of difference
      err_est = 2**mode / (2**mode - 1) * sqrt( sum( ssvd ** 2 ) )

      ! reset curret state
      call copy_basis(X%U(:rx), Utmp)
      X%S(:rx,:rx) = Stmp

      return
   end subroutine compute_splitting_error

   subroutine compute_increment_norm(nrm, U, S, Ulag, Slag, ifnorm)
      !! This function computes the norm of the solution increment in a cheap way avoiding the
      !! construction of the full low-rank solutions.
      real(wp),                               intent(out)   :: nrm
      !! Increment norm of current timestep
      class(abstract_vector_rdp),             intent(in)    :: U(:)
      !! Low-rank basis of current solution
      real(wp),                               intent(in)    :: S(:,:)
      !! Coefficients of current solution
      class(abstract_vector_rdp),             intent(in)    :: Ulag(:)
      !! Low-rank basis of lagged solution
      real(wp),                               intent(in)    :: Slag(:,:)
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
      call project_onto_common_basis(V1, V2, Ulag, U)

      ! project second low-rank state onto common basis and construct difference
      allocate(D(r+rl,r+rl)); D = 0.0_wp
      D(    :rl  ,     :rl  ) = Slag - matmul(V1, matmul(S, transpose(V1)))
      D(rl+1:rl+r,     :rl  ) =      - matmul(V2, matmul(S, transpose(V1)))
      D(    :rl  , rl+1:rl+r) =      - matmul(V1, matmul(S, transpose(V2)))
      D(rl+1:rl+r, rl+1:rl+r) =      - matmul(V2, matmul(S, transpose(V2)))

      ! compute Frobenius norm of difference
      nrm = sqrt(sum(svdvals(D))**2)
      if (ifnorm_) nrm = nrm/U(1)%get_size()

      return
   end subroutine compute_increment_norm

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
      real(wp),                               allocatable   :: wrk(:,:)
      integer                                               :: i, rk

      real(wp), allocatable :: ssvd_v(:)
      character(len=128) :: fmt
      integer :: j
      
      ! Optional argument
      trans = optval(iftrans, .false.)

      rk = X%rk
      allocate(R(rk,rk));   R = 0.0_wp 
      allocate(perm(rk));     perm = 0
      allocate(wrk(rk,rk)); wrk = 0.0_wp

      print *, ''
      ssvd_v = svdvals(X%S(:rk,:rk))
      print *, 'M'
      !write(fmt,'(A,I0,A)') '("SVD_M0  : ",', X%rk, '(F10.6,1X))'
      !print fmt, ssvd_v
      print *, ssvd_v(:4)

      ! Apply propagator to initial basis
      allocate(exptAU, source=X%U(1)); call exptAU%zero()
      do i = 1, rk
         call exptA(exptAU, A, X%U(i), tau, info, trans)
         call X%U(i)%axpby(0.0_wp, exptAU, 1.0_wp) ! overwrite old solution
      end do
      ! Reorthonormalize in-place
      call qr(X%U(:rk), R, perm, info)
      call check_info(info, 'qr_pivot', module=this_module, procedure='M_forward_map_rdp')

      call innerprod(wrk, X%U(:rk), X%U(:rk))
      wrk = wrk - eye(4)
      do i = 1, rk
         print '(4(E19.12,1X))', (wrk(i,j), j = 1, rk)
      end do

      ! Update low-rank fcators
      call apply_inverse_permutation_matrix(R, perm)

      print *,'R'
      ssvd_v = svdvals(R)
      print *, ssvd_v(:4)
      do i = 1, rk
         print '(4(E19.12,1X))', (R(i,j), j = 1, rk)
      end do
      print *,'S'
      ssvd_v = svdvals(X%S(:rk,:rk))
      print *, ssvd_v(:4)
      do i = 1, rk
         print '(4(E19.12,1X))', (X%S(i,j), j = 1, rk)
      end do

      ! Update coefficient matrix
      wrk = matmul(X%S(:rk,:rk), transpose(R))
      X%S(:rk,:rk) = matmul(R, wrk)

      ssvd_v = svdvals(X%S(:rk,:rk))
      !ssvd_v = svdvals(R)
      !write(fmt,'(A,I0,A)') '("SVD_R   : ",', X%rk, '(F10.6,1X))'
      !print fmt, ssvd_v
      print *, ssvd_v(:4)

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
      integer                                               :: rk, rkmax

      class(abstract_vector_rdp),  allocatable   :: U1(:)
      class(abstract_vector_rdp),  allocatable   :: BBTU(:)

      real(wp), allocatable :: ssvd_v(:)
      character(len=128) :: fmt
    
      rk = X%rk
      rkmax = size(X%U)
      if (.not. allocated(U1))   allocate(U1(  rkmax), source=X%U(1))
      if (.not. allocated(BBTU)) allocate(BBTU(rkmax), source=X%U(1))
      call zero_basis(U1); call zero_basis(BBTU)

      call K_step_lyapunov(X, U1(:rk), BBTU(:rk), B, tau, info)
      call S_step_lyapunov(X, U1(:rk), BBTU(:rk),    tau, info)
      call L_step_lyapunov(X, U1(:rk),            B, tau, info)
      
      ! Copy updated low-rank factors to output
      call copy_basis(X%U(:rk), U1(:rk))

      print *, 'G'
      ssvd_v = svdvals(X%S(:X%rk,:X%rk))
      !write(fmt,'(A,I0,A)') '("SVD_G   : ",', X%rk, '(F10.6,1X))'
      !print fmt, ssvd_v
      print *, ssvd_v(:4)
      print *, ''
               
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
      real(wp),                    allocatable   :: Swrk(:,:)

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
      real(wp),                    allocatable   :: Swrk(:,:)

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
      class(abstract_vector_rdp),  allocatable   :: Uwrk(:)

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

end module LightROM_LyapunovSolvers
