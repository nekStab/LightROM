module LightROM_SolverUtils
   ! stdlib
   use stdlib_optval, only : optval
   ! LightKrylov for Linear Algebra
   use LightKrylov
   use LightKrylov, only : dp
   use LightKrylov_Logger, only: log_message, log_information, log_warning, log_debug, check_info, stop_error
   ! LightROM
   use LightROM_AbstractLTIsystems
   use LightROM_LyapunovUtils
   use LightROM_RiccatiUtils
   use LightROM_Timing, only: lr_timer => global_lightROM_timer, time_lightROM
   
   implicit none 
   character(len=*), parameter, private :: this_module = 'LR_SolverUtils'

   ! global scratch arrays
   class(abstract_vector_rdp),  allocatable   :: Uwrk0(:)
   class(abstract_vector_rdp),  allocatable   :: Uwrk1(:)
   class(abstract_vector_rdp),  allocatable   :: U1(:)
   class(abstract_vector_rdp),  allocatable   :: QU(:)
   real(dp),                    allocatable   :: Swrk0(:,:)
   real(dp),                    allocatable   :: Swrk1(:,:)
   ! global scratch arrays for the predictor step
   !class(abstract_vector_rdp),  allocatable   :: U0(:)
   !class(abstract_vector_rdp),  allocatable   :: T0(:)
   !class(abstract_vector_rdp),  allocatable   :: Ut(:)
   !class(abstract_vector_rdp),  allocatable   :: Tt(:)
   !real(dp),                    allocatable   :: S0(:,:)
   ! global scratch arrays SVD
   real(dp),                    allocatable   :: ssvd(:)
   real(dp),                    allocatable   :: Usvd(:,:), VTsvd(:,:)

   public :: abstract_dlra_stepper
   public :: set_initial_rank_lyapunov_rdp, set_initial_rank_riccati_rdp
   public :: rank_adaptive_projector_splitting_DLRA_step
   public :: set_initial_rank
   !public :: M_map, G_map

   interface rank_adaptive_projector_splitting_DLRA_step
      module procedure rank_adaptive_projector_splitting_DLRA_step_lyapunov_rdp
      module procedure rank_adaptive_projector_splitting_DLRA_step_riccati_rdp
   end interface

   interface fixed_rank_projector_splitting_DLRA_step
      module procedure fixed_rank_projector_splitting_DLRA_step_lyapunov_rdp
      module procedure fixed_rank_projector_splitting_DLRA_step_riccati_rdp
   end interface

   interface set_initial_rank
      module procedure set_initial_rank_lyapunov_rdp
      module procedure set_initial_rank_riccati_rdp
   end interface

   interface M_map
      module procedure M_map_rdp
   end interface

   interface G_map
      module procedure G_map_lyapunov_rdp
      module procedure G_map_riccati_rdp
   end interface

   interface K_step
      module procedure K_step_lyapunov_rdp
      module procedure K_step_riccati_rdp
   end interface

   interface S_step
      module procedure S_step_lyapunov_rdp
      module procedure S_step_riccati_rdp
   end interface

   interface L_step
      module procedure L_step_lyapunov_rdp
      module procedure L_step_riccati_rdp
   end interface

   abstract interface
      subroutine abstract_dlra_stepper(if_rank_adaptive, info)
         logical, intent(in) :: if_rank_adaptive
         integer, intent(out) :: info
      end subroutine
   end interface

contains

   !---------------------------------------------------------------
   !
   !  Fixed-rank projector-splitting DLRA step
   !
   !---------------------------------------------------------------

   subroutine fixed_rank_projector_splitting_DLRA_step_lyapunov_rdp( &
            &  X, A, B, tau, mode, info, exptA, trans)
      !! Driver for the time-stepper defining the splitting logic for each step of the the 
      !! projector-splitting integrator
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_linop_rdp),              intent(inout) :: A
      !! Linear operator
      class(abstract_vector_rdp),             intent(in)    :: B(:)
      !! Low-Rank inhomogeneity.
      real(dp),                               intent(in)    :: tau
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
      character(len=*), parameter :: this_procedure = 'fixed_rank_projector_splitting_DLRA_step_lyapunov_rdp'
      character(len=128) :: msg
      integer :: istep, nsteps

      if (time_lightROM()) call lr_timer%start(this_procedure)

      select case (mode)
      case (1)
         ! Lie-Trotter splitting
         call M_map(X, A,     tau, info, exptA, trans)
         call G_map(X, B,     tau, info)
      case (2) 
         ! Strang splitting
         call M_map(X, A, 0.5*tau, info, exptA, trans)
         call G_map(X, B,     tau, info)
         call M_map(X, A, 0.5*tau, info, exptA, trans)
      end select

      if (time_lightROM()) call lr_timer%stop(this_procedure)
   end subroutine fixed_rank_projector_splitting_DLRA_step_lyapunov_rdp

   subroutine fixed_rank_projector_splitting_DLRA_step_riccati_rdp( &
            &  X, A, B, CT, Qc, Rinv, tau, mode, info, exptA, trans)
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
      real(dp),                               intent(in)    :: Qc(:,:)
      !! Measurement weights.
      real(dp),                               intent(in)    :: Rinv(:,:)
      !! Inverse of the actuation weights.
      real(dp),                               intent(in)    :: tau
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
      character(len=*), parameter :: this_procedure = 'fixed_rank_projector_splitting_DLRA_step_riccati_rdp'
      integer :: istep, nsteps, rk

      if (time_lightROM()) call lr_timer%start(this_procedure)

      select case (mode)
      case (1) 
         ! Lie-Trotter splitting
         call M_map(X, A,               tau, info, exptA, trans)
         call G_map(X, B, CT, Qc, Rinv, tau, info)
      case (2) 
         call stop_error('Second order integration not implemented for the Riccati equation', this_module, this_procedure)
         ! Strang splitting
         ! Prepare arrays for predictor step
         !rk = size(X%U)
         !! precomputation arrays
         !if (.not. allocated(U0)) allocate(U0(1:rk), source=X%U(1:rk))
         !if (.not. allocated(T0)) allocate(T0(1:rk), source=X%U(1:rk))
         !if (.not. allocated(Ut)) allocate(Ut(1:rk), source=X%U(1:rk))
         !if (.not. allocated(Tt)) allocate(Tt(1:rk), source=X%U(1:rk))
         !if (.not. allocated(S0)) allocate(S0(1:rk,1:rk))
         !call zero_basis(U0); call zero_basis(T0); call zero_basis(Ut); call zero_basis(Tt); S0 = 0.0_dp
         !! scratch arrays
         !if (.not. allocated(Uwrk0)) allocate(Uwrk0(1:rk), source=X%U(1:rk))
         !if (.not. allocated(Swrk0)) allocate(Swrk0(1:rk,1:rk))
         !call zero_basis(Uwrk0); Swrk0 = 0.0_dp
         !! second order step
         !call M_map        (X, A,                  0.5*tau, info, exptA, trans)
         !! Save current state
         !call copy(U0, X%U); S0 = X%S               ! --> save  
         !! Precompute T0
         !block
         !   class(abstract_vector_rdp), allocatable :: Xwrk(:)
         !   call linear_combination(Xwrk, U0, S0); call copy(Uwrk0, Xwrk) ! K0 = U0 @ S0
         !   call apply_premult_outerprod_w(Swrk0, U0, Uwrk0, B, Rinv)   ! (U0.T) @ B @ R^(-1) @ B.T @ K0
         !   call linear_combination(Xwrk, Uwrk0, Swrk0); call copy(T0, Xwrk) ! K0 @ Swrk0
         !end block
         !! First order integration
         !call G_map(X,    B, CT, Qc, Rinv,     tau, info, ifpred=.true., T0=T0)
         !! Precompute Tt
         !call copy(Ut, X%U)                                               ! --> save
         !block
         !   class(abstract_vector_rdp), allocatable :: Xwrk(:)
         !   call linear_combination(Xwrk, X%U, X%S); call copy(Uwrk0, Xwrk) ! Kt = Ut @ St
         !   call apply_premult_outerprod_w(Swrk0, X%U, Uwrk0, B, Rinv)  ! (Ut.T) @ B @ R^(-1) @ B.T @ Kt
         !   call linear_combination(Xwrk, Uwrk0, Swrk0); call copy(Tt, Xwrk) ! Kt @ Swrk0
         !end block
         !! Reset state
         !call copy(X%U, U0); X%S = S0
         !! Second order integration
         !call G_map(X,    B, CT, Qc, Rinv,     tau, info, ifpred=.false., &
         !                         & T0=T0, Tt=Tt, U0=U0, Ut=Ut)
         !call M_map        (X, A,                  0.5*tau, info, exptA, trans)
      end select

      if (time_lightROM()) call lr_timer%stop(this_procedure)
   end subroutine fixed_rank_projector_splitting_DLRA_step_riccati_rdp

   !---------------------------------------------------------------
   !
   !  Rank-adaptive projector-splitting DLRA step
   !
   !---------------------------------------------------------------

   subroutine generic_rank_adaptive_projector_splitting_DLRA_step_rdp( &
            &  X, info, rk_reduction_lock, tol, stepper, this_procedure)
      !! Wrapper for projector_splitting_DLRA_step_lyapunov_rdp adding the logic for rank-adaptivity
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      integer,                                intent(out)   :: info
      !! Information flag
      integer,                                intent(inout) :: rk_reduction_lock
      !! 'timer' to disable rank reduction
      real(dp),                               intent(in)    :: tol
      !! Tolerance for rank determination
      procedure(abstract_dlra_stepper)                      :: stepper
      !! Generic wrapper for the stepper (Lyapnuov/Riccati)
      character(len=*),                       intent(in)    :: this_procedure
      !! Specific procedure name (for error attribution)
      
      ! Internal variables
      integer                                               :: i, istep, rk, irk, rkmax, ndigits
      logical                                               :: accept_step, found
      real(dp),                               allocatable   :: coef(:)
      real(dp)                                              :: norm
      character(len=256)                                    :: msg, fmt
      integer, parameter                                    :: max_step = 40  ! might not be needed

      if (time_lightROM()) call lr_timer%start(this_procedure)

      ! ensure that we are integrating one more rank than we use for approximation
      X%rk = X%rk + 1
      rk = X%rk ! this is only to make the code more readable
      rkmax = size(X%U)
      ndigits = max(1,ceiling(log10(real(rkmax))))
      
      accept_step = .false.
      istep = 1
      do while ( .not. accept_step .and. istep < max_step )
         ! run a regular step
         call stepper(.false., info)
         if (info /= 0) then
            write(msg, '(A,I0,A)') 'Error in stepper at iteration ', i, '.'
            call stop_error(msg, this_module, this_procedure)
         end if
         ! compute singular values of X%S
         call svd(X%S(:rk,:rk), ssvd(:rk), Usvd(:rk,:rk), VTsvd(:rk,:rk))
         call find_rank(found, irk, ssvd(:rk), tol)
         
         ! choose action
         if (.not. found) then ! none of the singular values is below tolerance
            ! increase rank and run another step
            if (rk == rkmax) then ! cannot increase rank without reallocating X%U and X%S
               write(msg,'(A,I0,A,A)') 'Cannot increase rank, rkmax = ', rkmax, ' is reached. ', &
                        & 'Increase rkmax and restart!'
               call stop_error(msg, this_module, 'rank_adaptive_PS_DLRA_step_lyapunov_rdp')
            else
               write(fmt,'("(A,I3,A,I",I0,".",I0,",A,E14.8)")') ndigits, ndigits
               write(msg,fmt) 'rk= ', rk, ', s_', rk,' = ', ssvd(rk)
               call log_information(msg, this_module, 'DLRA_main')
               write(msg,'(A,I0)') 'Rank increased to rk= ', rk + 1
               call log_message(msg, this_module, 'DLRA_main')
               
               X%rk = X%rk + 1
               rk = X%rk ! this is only to make the code more readable
               ! set coefficients to zero (for redundancy)
               X%S(:rk, rk) = 0.0_dp 
               X%S( rk,:rk) = 0.0_dp
               ! add random vector ...
               call X%U(rk)%rand(.false.)
               ! ... and orthonormalize
               call orthogonalize_against_basis(X%U(rk), X%U(:rk-1), info, if_chk_orthonormal=.false.)
               call check_info(info, 'orthogonalize_against_basis', this_module, &
                                 & 'rank_adaptive_PS_DLRA_step_lyapunov_rdp')
               call X%U(rk)%scal(1.0_dp / X%U(rk)%norm())

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

      if (time_lightROM()) call lr_timer%stop(this_procedure)
   end subroutine generic_rank_adaptive_projector_splitting_DLRA_step_rdp

   subroutine rank_adaptive_projector_splitting_DLRA_step_lyapunov_rdp( &
            &  X, A, B, tau, mode, exptA, trans, info, rk_reduction_lock, tol)
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      class(abstract_linop_rdp),              intent(inout) :: A
      class(abstract_vector_rdp),             intent(in)    :: B(:)
      real(dp),                               intent(in)    :: tau
      integer,                                intent(in)    :: mode
      procedure(abstract_exptA_rdp)                         :: exptA
      logical,                                intent(in)    :: trans
      integer,                                intent(out)   :: info
      integer,                                intent(inout) :: rk_reduction_lock
      real(dp),                               intent(in)    :: tol
      ! internal
      character(len=*), parameter :: this_procedure = 'rank_adaptive_projector_splitting_DLRA_step_lyapunov_rdp'
      
      call generic_rank_adaptive_projector_splitting_DLRA_step_rdp( &
            &  X, info, rk_reduction_lock, tol, stepper, this_procedure)

   contains

      subroutine stepper(if_rank_adaptive, info)
         logical, intent(in) :: if_rank_adaptive
         integer, intent(out) :: info
         !! Information flag
         call fixed_rank_projector_splitting_DLRA_step_lyapunov_rdp( &
               &  X, A, B, tau, mode, info, exptA, trans)
      end subroutine

   end subroutine rank_adaptive_projector_splitting_DLRA_step_lyapunov_rdp

   subroutine rank_adaptive_projector_splitting_DLRA_step_riccati_rdp( &
            &  X, A, B, CT, Qc, Rinv, tau, mode, exptA, trans, info, rk_reduction_lock, tol)
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      class(abstract_linop_rdp),              intent(inout) :: A
      class(abstract_vector_rdp),             intent(in)    :: B(:)
      class(abstract_vector_rdp),             intent(in)    :: CT(:)
      real(dp),                               intent(in)    :: Qc(:,:)
      real(dp),                               intent(in)    :: Rinv(:,:)
      real(dp),                               intent(in)    :: tau
      integer,                                intent(in)    :: mode
      procedure(abstract_exptA_rdp)                         :: exptA
      logical,                                intent(in)    :: trans
      integer,                                intent(out)   :: info
      integer,                                intent(inout) :: rk_reduction_lock
      real(dp),                               intent(in)    :: tol
      ! internal
      character(len=*), parameter :: this_procedure = 'rank_adaptive_projector_splitting_DLRA_step_riccati_rdp'

      call generic_rank_adaptive_projector_splitting_DLRA_step_rdp( &
            &  X, info, rk_reduction_lock, tol, stepper, this_procedure)

   contains

      subroutine  stepper(if_rank_adaptive, info)
         logical, intent(in) :: if_rank_adaptive
         integer, intent(out) :: info
         !! Information flag
         call fixed_rank_projector_splitting_DLRA_step_riccati_rdp( &
               &  X, A, B, CT, Qc, Rinv, tau, mode, info, exptA, trans)
      end subroutine

   end subroutine rank_adaptive_projector_splitting_DLRA_step_riccati_rdp

   !---------------------------------------------------------------
   !
   !  Set initial rank
   !
   !---------------------------------------------------------------

   subroutine generic_set_initial_rank_rdp(X, tol, stepper, this_procedure, rk_init, nsteps)
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      real(dp),                               intent(in)    :: tol
      !! Tolerance on the last singular value to determine rank
      procedure(abstract_dlra_stepper)                      :: stepper
      !! Generic wrapper for the stepper (Lyapnuov/Riccati)
      character(len=*),                       intent(in)    :: this_procedure
      !! Specific procedure name (for error attribution)
      integer,                      optional, intent(in)    :: rk_init
      !! Smallest tested rank
      integer,                      optional, intent(in)    :: nsteps
      !! Number of steps to determine rank
      
      ! internal
      character(len=512) :: msg, fmt
      class(abstract_vector_rdp), allocatable :: Utmp(:)
      real(dp),                   allocatable :: Stmp(:,:), svals(:)
      integer :: i, n, rk, irk, info, rkmax
      logical :: found, accept_rank

      ! optional arguments
      X%rk = max(optval(rk_init, 1), 1)
      n = optval(nsteps, 5)
      rkmax = size(X%U)

      info = 0
      accept_rank = .false.

      ! save initial condition
      allocate(Utmp(rkmax), source=X%U)
      allocate(Stmp(rkmax,rkmax)); Stmp = X%S

      do while (.not. accept_rank .and. X%rk <= rkmax)
         write(msg,'(4X,A,I0)') 'Test r = ', X%rk
         call log_message(msg, this_module, this_procedure)
         ! run integrator
         do i = 1,n
            call stepper(.false., info)
            if (info /= 0) then
               write(msg, '(A,I0,A)') 'Error in stepper at iteration ', i, '.'
               call stop_error(msg, this_module, this_procedure)
            end if
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
         call log_debug(msg, this_module, this_procedure)
         if (found) then
            accept_rank = .true.
            X%rk = irk
            write(msg,'(4X,A,I2,A,E10.4)') 'Accpeted rank: r = ', X%rk-1, ',     s_{r+1} = ', svals(X%rk)
            call log_message(msg, this_module, this_procedure)
         else
            X%rk = 2*X%rk
         end if
         
         ! reset initial conditions
         call copy(X%U, Utmp)
         X%S = Stmp
      end do

      if (X%rk > rkmax) then
         write(msg, '(A)') 'Maximum rank reached but singular values are not converged. Increase rkmax and restart.'
         call stop_error(msg, this_module, this_procedure)
      end if

      ! reset to the rank of the approximation which we use outside of the integrator & mark rank as initialized
      X%rk = max(X%rk - 1, 1)
      X%rank_is_initialised = .true.
   end subroutine generic_set_initial_rank_rdp

   subroutine set_initial_rank_lyapunov_rdp(X, A, B, tau, mode, exptA, trans, tol, rk_init, nsteps)
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      class(abstract_linop_rdp),              intent(inout) :: A
      class(abstract_vector_rdp),             intent(in)    :: B(:)
      real(dp),                               intent(in)    :: tau
      integer,                                intent(in)    :: mode
      procedure(abstract_exptA_rdp)                         :: exptA
      logical,                                intent(in)    :: trans
      real(dp),                               intent(in)    :: tol
      integer,                      optional, intent(in)    :: rk_init
      integer,                      optional, intent(in)    :: nsteps
      ! internal
      character(len=*), parameter :: this_procedure = 'set_initial_rank_lyapunov_rdp'
      
      call generic_set_initial_rank_rdp(X, tol, stepper, this_procedure, rk_init, nsteps)

   contains

      subroutine stepper(if_rank_adaptive, info)
         logical, intent(in) :: if_rank_adaptive
         integer, intent(out) :: info
         call fixed_rank_projector_splitting_DLRA_step_lyapunov_rdp( &
               &  X, A, B, tau, mode, info, exptA, trans)
      end subroutine

   end subroutine set_initial_rank_lyapunov_rdp

   subroutine set_initial_rank_riccati_rdp(X, A, B, CT, Qc, Rinv, tau, mode, exptA, trans, tol, rk_init, nsteps)
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      class(abstract_linop_rdp),              intent(inout) :: A
      class(abstract_vector_rdp),             intent(in)    :: B(:)
      class(abstract_vector_rdp),             intent(in)    :: CT(:)
      real(dp),                               intent(in)    :: Qc(:,:)
      real(dp),                               intent(in)    :: Rinv(:,:)
      real(dp),                               intent(in)    :: tau
      integer,                                intent(in)    :: mode
      procedure(abstract_exptA_rdp)                         :: exptA
      logical,                                intent(in)    :: trans
      real(dp),                               intent(in)    :: tol
      integer,                      optional, intent(in)    :: rk_init
      integer,                      optional, intent(in)    :: nsteps
      ! internal
      character(len=*), parameter :: this_procedure = 'set_initial_rank_riccati_rdp'

      call generic_set_initial_rank_rdp(X, tol, stepper, this_procedure, rk_init, nsteps)

   contains

      subroutine stepper(if_rank_adaptive, info)
         logical, intent(in) :: if_rank_adaptive
         integer, intent(out) :: info
         call fixed_rank_projector_splitting_DLRA_step_riccati_rdp( &
               &  X, A, B, CT, Qc, Rinv, tau, mode, info, exptA, trans)
      end subroutine

   end subroutine set_initial_rank_riccati_rdp

   !---------------------------------------------------------------
   !
   !  M map
   !
   !---------------------------------------------------------------

   subroutine M_map_rdp(X, A, tau, info, exptA, iftrans)
      !! This subroutine computes the solution of the stiff linear part of the 
      !! differential equation exactly using the matrix exponential.
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_linop_rdp),              intent(inout) :: A
      !! Linear operator.
      real(dp),                               intent(in)    :: tau
      !! Time step.
      integer,                                intent(out)   :: info
      !! Information flag
      procedure(abstract_exptA_rdp)                         :: exptA
      !! Routine for computation of the exponential pabstract_vector),  ropagator (default: Krylov-based exponential operator).
      logical, optional,                      intent(in)    :: iftrans
      logical                                               :: trans
      !! Determine whether \(\mathbf{A}\) (default `.false.`) or \( \mathbf{A}^T\) (`.true.`) is used.

      ! Internal variables
      character(len=*), parameter :: this_procedure = 'M_map_rdp'
      class(abstract_vector_rdp),             allocatable   :: exptAU    ! scratch basis
      real(dp),                               allocatable   :: R(:,:)    ! QR coefficient matrix
      integer                                               :: i, rk

      if (time_lightROM()) call lr_timer%start(this_procedure)
      
      ! Optional argument
      trans = optval(iftrans, .false.)

      rk = X%rk
      allocate(R(rk,rk)); R = 0.0_dp

      ! Apply propagator to initial basis
      allocate(exptAU, source=X%U(1)); call exptAU%zero()
      do i = 1, rk
         call exptA(exptAU, A, X%U(i), tau, info, trans)
         call copy(X%U(i), exptAU) ! overwrite old solution
      end do
      ! Reorthonormalize in-place
      call qr(X%U(:rk), R, info)
      call check_info(info, 'qr', this_module, this_procedure)

      ! Update coefficient matrix
      X%S(:rk,:rk) = matmul(R, matmul(X%S(:rk,:rk), transpose(R)))

      if (time_lightROM()) call lr_timer%stop(this_procedure)
   end subroutine M_map_rdp

   !---------------------------------------------------------------
   !
   !  G map
   !
   !---------------------------------------------------------------

   subroutine G_map_lyapunov_rdp(X, B, tau, info)
      !! This subroutine computes the solution of the non-stiff part of the 
      !! differential equation numerically using first-order explicit Euler.
      !! The update of the full low-rank factorization requires three separate
      !! steps called K, S, L.
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_vector_rdp),             intent(in)    :: B(:)
      !! Low-Rank inhomogeneity.
      real(dp),                               intent(in)    :: tau
      !! Time step.
      integer,                                intent(out)   :: info
      !! Information flag.

      ! Internal variables
      character(len=*), parameter :: this_procedure = 'G_map_lyapunov_rdp'
      class(abstract_vector_rdp),  allocatable              :: U1(:)
      class(abstract_vector_rdp),  allocatable              :: BBTU(:)
      integer                                               :: rk

      if (time_lightROM()) call lr_timer%start(this_procedure)
    
      rk = X%rk
      allocate(  U1(rk), source=X%U(1)); call zero_basis(U1)
      allocate(BBTU(rk), source=X%U(1)); call zero_basis(BBTU)

      call K_step(X, U1, BBTU, B, tau, info)
      call S_step(X, U1, BBTU,    tau, info)
      call L_step(X, U1,       B, tau, info)
      
      ! Copy updated low-rank factors to output
      call copy(X%U(:rk), U1)

      if (time_lightROM()) call lr_timer%stop(this_procedure)
   end subroutine G_map_lyapunov_rdp

   subroutine G_map_riccati_rdp(X, B, CT, Qc, Rinv, tau, info, ifpred, T0, Tt, U0, Ut)
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
      real(dp),                               intent(in)    :: Qc(:,:)
      !! Measurement weights.
      real(dp),                               intent(in)    :: Rinv(:,:)
      !! Inverse of the actuation weights.
      real(dp),                               intent(in)    :: tau
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
      character(len=*), parameter :: this_procedure = 'G_map_riccati_rdp'
      integer :: rk

      if (time_lightROM()) call lr_timer%start(this_procedure)

      rk = size(X%U)
      if (.not. allocated(U1))  allocate(U1( 1:rk), source=X%U(1)); 
      if (.not. allocated(QU))  allocate(QU( 1:rk), source=X%U(1));
      if (.not. allocated(Swrk0)) allocate(Swrk0(1:rk,1:rk))
      call zero_basis(U1); call zero_basis(QU); Swrk0 = 0.0_dp

      if (present(ifpred)) then
         call stop_error('Second order integration not implemented for the Riccati equation', this_module, this_procedure)
         ! second order in time
         !if (ifpred) then ! predictor step with precomputed T0
         !   call K_step(X, U1, QU, B, CT, Qc, Rinv,     tau, info, reverse=.false., NL=T0)
         !   call S_step(X, U1, QU, B, CT, Qc, Rinv,     tau, info, reverse=.false.)
         !   call L_step(X, U1,     B, CT, Qc, Rinv,     tau, info)
         !else             ! corrector step with precomputed T0, Tt and U0, Ut
         !   ! forward steps based on T0, U0 (within X)
         !   call K_step(X, U1, QU, B, CT, Qc, Rinv, 0.5*tau, info, reverse=.false., NL=T0)
         !   call S_step(X, U1, QU, B, CT, Qc, Rinv, 0.5*tau, info, reverse=.false.)
         !   call L_step(X, U1,     B, CT, Qc, Rinv,     tau, info)
         !   ! Compute Gamma = 0.5*(T0 @ (U1.T @ U0) + Tt @ (U1.T @ Ut))
         !   call zero_basis(QU); Swrk0 = 0.0_dp                         ! we use QU as a scratch array
         !   block
         !      class(abstract_vector_rdp), allocatable :: Xwrk(:)
         !      Swrk0 = innerprod(X%U, U0)
         !      call linear_combination(Xwrk, T0, Swrk0); call copy(QU, Xwrk)
         !      Swrk0 = innerprod(X%U, Ut)
         !      call linear_combination(Xwrk, Tt, Swrk0); call copy(T0, Xwrk) ! overwrite T0 with Gamma
         !   end block
         !   call axpby_basis(0.5_dp, QU, 0.5_dp, T0)
         !   ! Update X to most recent value
         !   call copy(X%U, U1)
         !   ! reverse steps based on Gamma
         !   call S_step(X, U1, QU, B, CT, Qc, Rinv, 0.5*tau, info, reverse=.true., NL=T0)
         !   call K_step(X, U1, QU, B, CT, Qc, Rinv, 0.5*tau, info, reverse=.true., NL=T0)
         !end if
      else
         ! first order in time
         call K_step(   X, U1, QU, B, CT, Qc, Rinv,     tau, info)
         call S_step(   X, U1, QU, B, CT, Qc, Rinv,     tau, info)
         call L_step(   X, U1,     B, CT, Qc, Rinv,     tau, info)
      end if
      
      ! Copy updated low-rank factors to output
      call copy(X%U, U1)

      if (time_lightROM()) call lr_timer%stop(this_procedure)
   end subroutine G_map_riccati_rdp

   !---------------------------------------------------------------
   !
   !  K step
   !
   !---------------------------------------------------------------

   subroutine K_step_lyapunov_rdp(X, U1, BBTU, B, tau, info)
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_vector_rdp),             intent(out)   :: U1(:)
      !! Intermediate low-rank factor.
      class(abstract_vector_rdp),             intent(out)   :: BBTU(:)
      !! Precomputed application of the inhomogeneity.
      class(abstract_vector_rdp),             intent(in)    :: B(:)
      !! Low-Rank inhomogeneity.
      real(dp),                               intent(in)    :: tau
      !! Time step.
      integer,                                intent(out)   :: info
      !! Information flag.

      ! Internal variables
      character(len=*), parameter :: this_procedure = 'K_step_lyapunov_rdp'
      class(abstract_vector_rdp), allocatable :: Uwrk(:)
      integer :: rk

      if (time_lightROM()) call lr_timer%start(this_procedure)

      rk = X%rk
      call linear_combination(Uwrk, X%U(:rk), X%S(:rk,:rk))  ! K0
      call copy(U1, Uwrk)
      call apply_outerprod(BBTU, B, X%U(:rk))                ! Kdot
      ! Construct intermediate solution U1
      call axpby_basis(tau, BBTU, 1.0_dp, U1)                ! K0 + tau*Kdot
      ! Orthonormalize in-place
      call qr(U1, X%S(:rk,:rk), info)
      call check_info(info, 'qr', this_module, this_procedure)

      if (time_lightROM()) call lr_timer%stop(this_procedure)
   end subroutine K_step_lyapunov_rdp

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
      real(dp),                              intent(in)     :: Qc(:,:)
      !! Measurement weights.
      real(dp),                              intent(in)     :: Rinv(:,:)
      !! Inverse of the actuation weights.
      real(dp),                              intent(in)     :: tau
      !! Time step.
      integer,                               intent(out)    :: info
      !! Information flag.
      logical, optional,                     intent(in)     :: reverse
      !! For Strang splitting: Determine if we are in forward or reverse branch
      class(abstract_vector_rdp), optional,  intent(in)     :: NL(:)
      !! Precomputed non-linear term.

      ! Internal variables
      character(len=*), parameter :: this_procedure = 'K_step_riccati_rdp'
      logical :: reverse_order
      integer :: rk

      if (time_lightROM()) call lr_timer%start(this_procedure)

      ! Optional arguments
      reverse_order = optval(reverse, .false.)

      info = 0

      rk = size(X%U)
      if (.not. allocated(Uwrk0)) allocate(Uwrk0(1:rk), source=X%U(1)); 
      if (.not. allocated(Swrk0)) allocate(Swrk0(1:rk,1:rk)); 
      call zero_basis(Uwrk0); Swrk0 = 0.0_dp

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
      call axpby_basis(1.0_dp, QU, -1.0_dp, Uwrk0)

      ! Construct intermediate solution U1
      call axpby_basis(tau, Uwrk0, 1.0_dp, U1)                   ! K0 + tau*Kdot

      ! Orthonormalize in-place
      call qr(U1, Swrk0, info)
      call check_info(info, 'qr', this_module, this_procedure)
      X%S = Swrk0

      if (time_lightROM()) call lr_timer%stop(this_procedure)
   end subroutine K_step_riccati_rdp

   !---------------------------------------------------------------
   !
   !  S step
   !
   !---------------------------------------------------------------

   subroutine S_step_lyapunov_rdp(X, U1, BBTU, tau, info)
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_vector_rdp),             intent(in)    :: U1(:)
      !! Intermediate low-rank factor.
      class(abstract_vector_rdp),             intent(in)    :: BBTU(:)
      !! Precomputed application of the inhomogeneity.
      real(dp),                               intent(in)    :: tau
      !! Time step.
      integer,                                intent(out)   :: info
      !! Information flag.

      ! Internal variables
      character(len=*), parameter :: this_procedure = 'S_step_lyapunov_rdp'
      real(dp), allocatable :: Swrk(:,:)
      integer :: rk

      if (time_lightROM()) call lr_timer%start(this_procedure)

      rk = X%rk
      allocate(Swrk(rk,rk)); Swrk = 0.0_dp
      Swrk = innerprod(U1, BBTU)          ! - Sdot
      ! Construct intermediate coefficient matrix
      X%S(:rk,:rk) = X%S(:rk,:rk) - tau*Swrk

      if (time_lightROM()) call lr_timer%stop(this_procedure)
   end subroutine S_step_lyapunov_rdp

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
      real(dp),                              intent(in)    :: Qc(:,:)
      !! Measurement weights.
      real(dp),                              intent(in)    :: Rinv(:,:)
      !! Inverse of the actuation weights.
      real(dp),                              intent(in)    :: tau
      !! Time step.
      integer,                               intent(out)   :: info
      !! Information flag.
      logical, optional,                     intent(in)    :: reverse
      !! For Strang splitting: Determine if we are in forward or reverse branch
      class(abstract_vector_rdp), optional,  intent(in)    :: NL(:)
      !! Precomputed non-linear term.

      ! Internal variables
      character(len=*), parameter :: this_procedure = 'S_step_riccati_rdp'
      logical :: reverse_order
      integer :: rk

      if (time_lightROM()) call lr_timer%start(this_procedure)

      info = 0

      ! Optional arguments
      reverse_order = optval(reverse, .false.)

      rk = size(X%U)
      if (.not. allocated(Swrk0)) allocate(Swrk0(1:rk,1:rk))
      if (.not. allocated(Swrk1)) allocate(Swrk1(1:rk,1:rk))
      Swrk0 = 0.0_dp; Swrk1 = 0.0_dp

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

      if (time_lightROM()) call lr_timer%stop(this_procedure)
   end subroutine S_step_riccati_rdp

   !---------------------------------------------------------------
   !
   !  L step
   !
   !---------------------------------------------------------------

   subroutine L_step_lyapunov_rdp(X, U1, B, tau, info)
      class(abstract_sym_low_rank_state_rdp), intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_vector_rdp),             intent(in)    :: U1(:)
      !! Intermediate low-rank factor (from K step).
      class(abstract_vector_rdp),             intent(in)    :: B(:)
      !! Low-Rank inhomogeneity.
      real(dp),                               intent(in)    :: tau
      !! Time step.
      integer,                                intent(out)   :: info
      !! Information flag.

      ! Internal variables
      character(len=*), parameter :: this_procedure = 'L_step_lyapunov_rdp'
      class(abstract_vector_rdp), allocatable :: Uwrk(:)
      integer :: rk

      if (time_lightROM()) call lr_timer%start(this_procedure)

      rk = X%rk
      call linear_combination(Uwrk, X%U(:rk), transpose(X%S(:rk,:rk)))  ! L0.T
      ! Construct derivative
      call apply_outerprod(X%U(:rk), B, U1)       ! Ldot.T
      ! Construct solution L1.T
      call axpby_basis(tau, X%U(:rk), 1.0_dp, Uwrk)
      ! Update coefficient matrix
      X%S(:rk,:rk) = innerprod(Uwrk, U1)

      if (time_lightROM()) call lr_timer%stop(this_procedure)
   end subroutine L_step_lyapunov_rdp

   subroutine L_step_riccati_rdp(X, U1, B, CT, Qc, Rinv, tau, info)
      class(abstract_sym_low_rank_state_rdp),intent(inout) :: X
      !! Low-Rank factors of the solution.
      class(abstract_vector_rdp),            intent(in)    :: U1(:)
      !! Intermediate low-rank factor.
      class(abstract_vector_rdp),            intent(in)    :: B(:)
      !! System input.
      class(abstract_vector_rdp),            intent(in)    :: CT(:)
      !! System output.
      real(dp),                              intent(in)    :: Qc(:,:)
      !! Measurement weights.
      real(dp),                              intent(in)    :: Rinv(:,:)
      !! Inverse of the actuation weights.
      real(dp),                              intent(in)    :: tau
      !! Time step.
      integer,                               intent(out)   :: info
      !! Information flag.

      ! Internal variables
      character(len=*), parameter :: this_procedure = 'L_step_riccati_rdp'
      integer :: rk

      if (time_lightROM()) call lr_timer%start(this_procedure)

      info = 0

      rk = size(X%U)
      if (.not. allocated(Uwrk0)) allocate(Uwrk0(1:rk), source=X%U(1))
      if (.not. allocated(Uwrk1)) allocate(Uwrk1(1:rk), source=X%U(1))
      if (.not. allocated(Swrk0)) allocate(Swrk0(1:rk,1:rk)) 
      call zero_basis(Uwrk0); call zero_basis(Uwrk1); Swrk0 = 0.0_dp

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
      call axpby_basis(-1.0_dp, X%U, 1.0_dp, Uwrk0)

      ! Construct solution L1.T
      call axpby_basis(tau, Uwrk0, 1.0_dp, Uwrk1)               ! L0.T + tau*Ldot.T

      ! Update coefficient matrix
      X%S = innerprod(Uwrk1, U1)

      if (time_lightROM()) call lr_timer%stop(this_procedure)
   end subroutine L_step_riccati_rdp

end module LightROM_SolverUtils
