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

   ! module name
   private :: this_module
   character*128, parameter :: this_module = 'LightKrylov_LyapunovSolvers'
   ! rank_reduction_lock
   integer, private :: rk_reduction_lock

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
                                                                    & exptA, iftrans, ifverb, ifrk, tol)
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
      real(wp),                      optional, intent(in)   :: tol
      real(wp)                                              :: tolerance
      !! tolerance for the smallest singular value (this value is only used if ifrk = .true.)

      ! Internal variables   
      integer                                               :: istep, nsteps, iostep
      real(wp)                                              :: T
      
      integer                                               :: mode
      procedure(abstract_exptA_rdp), pointer                :: p_exptA => null()

      ! Optional arguments
      trans        = optval(iftrans, .false.)
      verbose      = optval(ifverb, .false.)
      rankadaptive = optval(ifrk, .false.)
      tolerance    = optval(tol, 1e-8_wp)
      if (present(exptA)) then
         p_exptA => exptA
      else
         p_exptA => k_exptA_rdp
      end if

      ! Initialize
      T                 = 0.0_wp
      rk_reduction_lock = 0

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
         write(*,*) "INFO : Time-integration order for the operator splitting of d > 2 &
                     &requires adjoint solves and is not implemented."
         write(*,*) "       Resetting torder = 2." 
         info = 1
         mode = 2
      else 
         write(*,*) "INFO : Invalid time-integration order specified."
         info = -1
         mode = 2
      endif

      dlra : do istep = 1, nsteps
         ! dynamical low-rank approximation solver
         if (rankadaptive) then
            call DLRA_rank_adaptive_lyapunov_step_rdp(X, A, B, tau, mode, info, p_exptA, tol, trans)
         else
            call numerical_low_rank_splitting_lyapunov_step_rdp(X, A, B, tau, mode, info, p_exptA, trans)
         end if

         T = T + tau
         ! here we can do some checks such as whether we have reached steady state
         if ( mod(istep,iostep) .eq. 0 ) then
            if (verbose) then
               write(*, *) "INFO : ", ISTEP, " steps of DLRA computed. T = ",T
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

   subroutine numerical_low_rank_splitting_lyapunov_step_rdp(X, A, B, tau, mode, info, exptA, iftrans)
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
      logical,                      optional, intent(in)    :: iftrans
      !! Determine whether \(\mathbf{A}\) (default `.false.`) or \( \mathbf{A}^T\) (`.true.`) is used.
      
      ! Internal variables
      integer                                               :: istep, nsteps
      logical                                               :: trans

      ! Optional argument
      trans = optval(iftrans, .false.)

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
      allocate(R(1:rk,1:rk));   R = 0.0_wp 
      allocate(perm(1:rk));     perm = 0
      allocate(wrk(1:rk,1:rk)); wrk = 0.0_wp

      ! Apply propagator to initial basis
      allocate(exptAU, source=X%U(1)); call exptAU%zero()
      do i = 1, rk
         call exptA(exptAU, A, X%U(i), tau, info, trans)
         call X%U(i)%axpby(0.0_wp, exptAU, 1.0_wp) ! overwrite old solution
      end do
      ! Reorthonormalize in-place
      call qr(X%U(1:rk), R, perm, info)
      call check_info(info, 'qr_pivot', module=this_module, procedure='M_forward_map_rdp')
      ! Update low-rank fcators
      call apply_inverse_permutation_matrix(R, perm)
      ! Update coefficient matrix
      wrk = matmul(X%S(1:rk,1:rk), transpose(R))
      X%S(1:rk,1:rk) = matmul(R, wrk)

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
      if (.not. allocated(U1))   allocate(U1(  1:rkmax), source=X%U(1))
      if (.not. allocated(BBTU)) allocate(BBTU(1:rkmax), source=X%U(1))
      call zero_basis(U1); call zero_basis(BBTU)

      call K_step_lyapunov(X, U1(1:rk), BBTU(1:rk), B, tau, info)
      call S_step_lyapunov(X, U1(1:rk), BBTU(1:rk),    tau, info)
      call L_step_lyapunov(X, U1(1:rk),             B, tau, info)
      
      ! Copy updated low-rank factors to output
      call copy_basis(X%U(1:rk), U1(1:rk))
               
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
      if (.not. allocated(Swrk)) allocate(Swrk(1:rkmax,1:rkmax)); Swrk = 0.0_wp
      allocate(perm(1:rk)); perm = 0
      block
         class(abstract_vector_rdp), allocatable :: Xwrk(:)
         call linear_combination(Xwrk, X%U(1:rk), X%S(1:rk,1:rk))             ! K0
         call copy_basis(U1, Xwrk)
      end block
      call apply_outerprod(BBTU, B, X%U(1:rk))   ! Kdot
      ! Construct intermediate solution U1
      call axpby_basis(U1, 1.0_wp, BBTU(1:rk), tau)   ! K0 + tau*Kdot
      ! Orthonormalize in-place
      call qr(U1(1:rk), Swrk(1:rk,1:rk), perm, info)
      call check_info(info, 'qr_pivot', module=this_module, procedure='K_step_Lyapunov_rdp')
      call apply_inverse_permutation_matrix(Swrk(1:rk,1:rk), perm)
      X%S(1:rk,1:rk) = Swrk(1:rk,1:rk)

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
      if (.not. allocated(Swrk)) allocate(Swrk(1:rkmax,1:rkmax)); Swrk = 0.0_wp
      call innerprod(Swrk(1:rk,1:rk), U1, BBTU)          ! - Sdot
      ! Construct intermediate coefficient matrix
      X%S(1:rk,1:rk) = X%S(1:rk,1:rk) - tau*Swrk(1:rk,1:rk)

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
      if (.not. allocated(Uwrk)) allocate(Uwrk(1:rkmax), source=X%U(1)); call zero_basis(Uwrk)

      block
         class(abstract_vector_rdp), allocatable :: Xwrk(:)
         call linear_combination(Xwrk, X%U(1:rk), transpose(X%S(1:rk,1:rk)))  ! L0.T
         call copy_basis(Uwrk(1:rk), Xwrk)
      end block
      call apply_outerprod(X%U(1:rk), B, U1)       ! Ldot.T
      ! Construct solution L1.T
      call axpby_basis(Uwrk(1:rk), 1.0_wp, X%U(1:rk), tau)
      ! Update coefficient matrix
      call innerprod(X%S(1:rk,1:rk), Uwrk(1:rk), U1)

      return
   end subroutine L_step_lyapunov_rdp

   !-----------------------------
   !
   !     RANK-ADAPTIVE PSI 
   !
   !-----------------------------

   subroutine DLRA_rank_adaptive_lyapunov_step_rdp(X, A, B, tau, mode, info, exptA, tol, iftrans)
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
      real(wp),                               intent(in)    :: tol
      !! tolerance for the smallest singular value
      logical,                      optional, intent(in)    :: iftrans
      logical                                               :: trans
      !! Determine whether \(\mathbf{A}\) (default `.false.`) or \( \mathbf{A}^T\) (`.true.`) is used.
      
      ! Internal variables
      integer                                               :: istep, irk
      logical                                               :: accept_step, found
      real(wp),                               allocatable   :: coef(:,:)
      real(wp)                                              :: norm
      ! svd
      real(wp),                               allocatable   :: svals(:)
      real(wp),                               allocatable   :: U(:,:), VT(:,:)
      character*128                                         :: msg

      integer, parameter                                    :: max_step = 5  ! might not be needed

      ! Optional argument
      trans = optval(iftrans, .false.)

      ! run a regular step with rk + 1 ranks
      X%rk = X%rk + 1
      call numerical_low_rank_splitting_lyapunov_step_rdp(X, A, B, tau, mode, info, exptA, trans)

      accept_step = .false.
      istep = 1
      do while ( .not. accept_step .and. istep < max_step )
         ! compute singular values of X%S
         allocate(U(1:X%rk,1:X%rk)); allocate(svals(1:X%rk)); allocate(VT(1:X%rk,1:X%rk)); 
         U = 0.0_wp; svals = 0.0_wp; VT = 0.0_wp
         call svd(X%S(1:X%rk, 1:X%rk), U, svals, VT)
         ! find singular value that falls below tolerance (if any)
         found = .false.
         tol_chk: do irk = 1, X%rk
            if ( svals(irk) < tol ) then
               found = .true.
               exit tol_chk
            end if
         end do tol_chk
         
         ! choose action
         if (.not. found) then ! none of the singular values is below tolerance
            ! increase rank and run another step
            if (X%rk == size(X%U)) then ! cannot increase rank without reallocating X%U and X%S
               info = -2
               write(msg, *) 'Tried to increase rank but rkmax is reached. Increase rkmax and restart!'
               call logger%log_error('Rank increase', module=this_module, procedure='DLRA_rank_adaptive_lyapunov_step_rdp', &
                                    & stat=info, errmsg=msg)
            else
               
               X%rk = X%rk + 1
               ! set coefficients to zero (for redundancy)
               X%S(1:X%rk, X%rk) = 0.0_wp 
               X%S(X%rk, 1:X%rk) = 0.0_wp
               ! add random vector ...
               call X%U(X%rk)%rand(.false.)
               ! ... and orthonormalize
               allocate(coef(X%rk,1)); coef = 0.0_wp
               call innerprod(coef, X%U(1:X%rk-1), X%U(X%rk:X%rk))
               block
                  class(abstract_vector_rdp), allocatable :: proj(:)
                  call linear_combination(proj, X%U(1:X%rk), coef)
                  call X%U(X%rk)%sub(proj(1))
               end block
               call X%U(X%rk)%scal(1.0_wp / X%U(X%rk)%norm())

               rk_reduction_lock = 10 ! avoid rank oscillations

            end if
         else ! the rank of the solution is sufficient

            accept_step = .true.

            if (irk /= X%rk .and. rk_reduction_lock == 0) then ! we should decrease the rank
               ! decrease rank

               ! rotate basis onto principal axes
               block
                  class(abstract_vector_rdp), allocatable :: Xwrk(:)
                  call linear_combination(Xwrk, X%U(1:X%rk), U)
                  call copy_basis(X%U(1:X%rk), Xwrk)
               end block
               X%S(1:X%rk, 1:X%rk) = diag(svals(1:X%rk))

               X%rk = max(irk, X%rk - 2)  ! reduce by at most 2              

            end if

         end if ! found
         istep = istep + 1
      end do ! while .not. accept_step

      ! decrease rk_reduction_lock
      if (rk_reduction_lock > 0) rk_reduction_lock = rk_reduction_lock - 1
      
      ! reset to the rank of the approximation which we use outside of the integrator 
      X%rk = X%rk - 1      

      return
   end subroutine DLRA_rank_adaptive_lyapunov_step_rdp

end module LightROM_LyapunovSolvers
