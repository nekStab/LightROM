module LightROM_RiccatiSolvers
   use LightKrylov
   use LightKrylov_expmlib
   use lightkrylov_BaseKrylov
   use LightROM_LyapunovUtils
   use LightROM_LyapunovSolvers
   use LightROM_utils
   use lightrom_AbstractLTIsystems
   use stdlib_optval, only : optval
   use stdlib_linalg, only : eye
   implicit none

   !> work arrays
   class(abstract_vector),  allocatable   :: Uwrk0(:)
   class(abstract_vector),  allocatable   :: Uwrk1(:)
   class(abstract_vector),  allocatable   :: QU(:)
   real(kind=wp),           allocatable   :: Swrk0(:,:)
   real(kind=wp),           allocatable   :: Swrk1(:,:)
   real(kind=wp),           allocatable   :: BTU(:,:)
   real(kind=wp),           allocatable   :: UTB(:,:)
   integer,                 allocatable   :: perm(:)

   private
   public :: numerical_low_rank_splitting_integrator_riccati
   public :: G_forward_map_riccati, K_step_riccati, S_step_riccati, L_step_riccati

   contains

   !------------------------------------------------
   !-----                                      -----
   !-----     Low-rank splitting integrator    -----
   !-----                                      -----
   !------------------------------------------------
   
   !=======================================================================================
   ! Numerical Low-Rank Splitting Integrator Subroutine
   !=======================================================================================
   !
   ! Purpose:
   ! --------
   ! Implementation of a Numerical Low-Rank Splitting Integrator as in Mena et al. (2018).
   !
   ! The algorithm constitutes a numerical integrator for the matrix-valued differential Riccati
   ! equation of the form
   !
   !          dX/dt = A @ X + X @ A.T + C.T @ Q @ C - X @ B @ R^{-1} @ B.T @ X .
   !
   ! where A is a nxn Hurwitz matrix, X is SPD and B and C.T are low-rank matrices.
   !
   ! Since A is Hurwitz, the equations converges to steady state for t -> infty, which corresponds
   ! to the associated algebrait Riccati equation of the form
   !
   !              0 = A @ X + X @ A.T + C.T @ Q @ C - X @ B @ R^{-1} @ B.T @ X .
   !
   ! The algorithm is based on three ideas:
   ! a) The operator splitting scheme proposed by Lubich & Oseledets (2014) that splits the 
   !    right-hand side of the differential equation into a linear stiff part that is solved
   !    explicitly and a possibly non-linear non-stiff part which is solved numerically. The
   !    two operators are then composed to obtain the integrator for the full Riccati equation.
   ! b) The Dynamic Low-Rank Approximation for the solution of general matrix differential 
   !    equations proposed by Nonnenmacher & Lubich (2007) which seeks to integrate only the
   !    leading low-rank factors of the solution to a large system by updating the matrix 
   !    factorization. The dynamical low-rank approximation scheme for the low-rank factors 
   !    of the solution is itself solved using a projector-splitting technique to cheaply 
   !    maintain orthonormality or the low-rank basis without explicit SVDs. 
   ! c) This algorithm has been applied to a large scale system by Mena et al. (2018) and we
   !    follow the algorithm proposed there.
   !
   ! Algorithmic Features:
   ! ----------------------
   ! - Separate integration of the stiff inhomogeneous part of the Riccati equation and the
   !   non-stiff inhomogeneity
   ! - Rank preserving time-integration that maintains orthonormality of the factorization
   !   basis
   ! - The stiff part of the problem is solved using a time-stepper approach to approximate 
   !   the action of the exponential propagator
   !
   ! Advantages:
   ! -----------
   ! - Rank of the approximate solution is user defined
   ! - The timesteps of the stiff and non-stiff parts of the code are independent
   ! - The integrator is adjoint-free
   ! - The operator of the homogeneous part and the inhomogeneity are not needed explicitly
   !   i.e. the algorithm is amenable to solution using Krylov methods (in particular for 
   !   the solution of the stiff part of the problem).
   ! - No SVDs are necessary for this alogorithm.
   !
   ! Limitations:
   ! ------------
   ! - Rank of the approximate solution is user defined. The appropriateness of this 
   !   approximation is not considered
   ! - The current implementation does not require an adjoint integrator. This means that
   !   the temporal order of the basic operator splitting scheme is limited to 1 (Lie-Trotter
   !   splitting) or at most 2 (Strang splitting). Higher order integrators are possible, but 
   !   require at least some backward integration (via the adjoint) in BOTH parts of the splitting. 
   !   (see Sheng-Suzuki and Goldman-Kaper theorems)
   !
   ! Input/Output Parameters:
   ! ------------------------
   ! - U      : Low-Rank factor of the solution           [Input/Output]
   ! - S      : Coefficients of the solution              [Input/Output]
   ! - A      : Linear Operator (SPD)                     [Input]
   ! - B      : Low-rank input for the Riccati equation   [Input]
   ! - CT : Low-rank input for the Riccati equation   [Input]
   ! - Tend   : Time horizon for integration              [Input]
   ! - tau    : Integration time-step                     [Input]
   ! - torder : Order of the time integration (1/2)       [Input]
   ! - info   : Iteration Information flag                [Output]
   ! - tol    : Tolerance for convergence                 [Optional, Input]
   ! - exptA  : Procedure for exponential propagator      [Optiomal, Input}
   !
   ! References:
   ! -----------
   ! Koch, O.,Lubich, C. (2007). "Dynamical Low‐Rank Approximation", SIAM Journal on Matrix Analysis 
   !     and Applications 29(2), 434-454
   ! Lubich, C., Oseledets, I.V. (2014). "A projector-splitting integrator for dynamical low-rank 
   !     approximation", BIT Numerical Mathematics 54, 171–188
   ! Mena, H., Ostermann, A., Pfurtscheller, L.-M., Piazzola, C. (2018). "Numerical low-rank 
   !     approximation of matrix differential equations", Journal of Computational and Applied Mathematics,
   !     340, 602-614
   !
   !=======================================================================================
   subroutine numerical_low_rank_splitting_integrator_riccati(U,S,LTI,Q,Rinv,Tend,tau,torder,info,exptA)
      !> Low-rank factors
      class(abstract_vector),              intent(inout) :: U(:) ! basis
      real(kind=wp),                       intent(inout) :: S(:,:) ! coefficients
      ! LTI system
      class(abstract_lti_system),          intent(in)    :: LTI
      class(abstract_linop),               intent(in)    :: Q
      real(kind=wp),                       intent(in)    :: Rinv(:,:)
      !> Integration time and step size
      real(kind=wp),                       intent(in)    :: Tend
      real(kind=wp),                       intent(inout) :: tau  ! desired/effective time step
      !> Order of time integration
      integer,                             intent(in)    :: torder
      !> Information flag
      integer,                             intent(out)   :: info
      !> Routine for computation of the exponential propagator
      procedure(abstract_exptA), optional                :: exptA

      !> Local variables   
      integer                                :: istep, nsteps, iostep
      real(kind=wp)                          :: T
      logical, parameter                     :: verbose = .false.
      procedure(abstract_exptA), pointer     :: p_exptA => null()

      T = 0.0_wp

      ! Optional arguments
      if (present(exptA)) then
         p_exptA => exptA
      else
         p_exptA => k_exptA
      endif

      ! --> Compute number of steps
      nsteps = floor(Tend/tau)

      iostep = nsteps/10
      if ( iostep .eq. 0 ) then
         iostep = 10
      endif

      dlra : do istep = 1, nsteps
         
         !> dynamical low-rank approximation step
         call numerical_low_rank_splitting_riccati_step(U, S, LTI, Q, Rinv, tau, torder, info, p_exptA)

         T = T + tau
         !> here we should do some checks such as whether we have reached steady state
         if ( mod(istep,iostep) .eq. 0 ) then
            if (verbose) then
               write(*, *) "INFO : ", ISTEP, " steps of DLRA computed. T = ",T
            endif
         endif
      enddo dlra

      deallocate(BTU); deallocate(UTB)

      return

   end subroutine numerical_low_rank_splitting_integrator_riccati

   !-----------------------------
   !-----     UTILITIES     -----
   !-----------------------------
   subroutine numerical_low_rank_splitting_riccati_step(U, S, LTI, Q, Rinv, tau, torder, info, exptA)
      !> Low-rank factors
      class(abstract_vector),     intent(inout) :: U(:) ! basis
      real(kind=wp),              intent(inout) :: S(:,:) ! coefficients
      ! LTI system
      class(abstract_lti_system), intent(in)    :: LTI
      class(abstract_linop),      intent(in)    :: Q
      real(kind=wp),              intent(in)    :: Rinv(:,:)
      !> Integration step size
      real(kind=wp),              intent(in)    :: tau
      !> Order of time integration
      integer,                    intent(in)    :: torder
      !> Information flag
      integer,                    intent(out)   :: info
      !> Routine for computation of the exponential propagator
      procedure(abstract_exptA)                 :: exptA

      !> Local variables
      integer                                   :: istep, nsteps, integrator, rk
      class(abstract_vector), allocatable       :: Uwrk0(:) ! basis
      real(kind=wp),          allocatable       :: Swrk0(:,:) ! coefficients

      if ( torder .eq. 1 ) then 
         integrator = 1
      else if ( torder .eq. 2 ) then
         integrator = 2
      else if ( torder .gt. 2 ) then
         write(*,*) "INFO : Time-integration order for the operator splitting of d > 2 &
                     &requires adjoint solves and is not implemented."
         write(*,*) "       Resetting torder = 2." 
         info = 1
         integrator = 2
      else 
         write(*,*) "INFO : Invalid time-integration order specified."
         info = -1
         integrator = 2
      endif

      select case (integrator)
      case (1) ! Lie-Trotter splitting
         call M_forward_map        (U,     S,     LTI,             tau, info, exptA, iftrans=.true.)
         call G_forward_map_riccati(U,     S,     LTI, Q, Rinv,    tau, info)
      case (2) ! Strang splitting
         rk = size(U)
         if (.not. allocated(Uwrk0)) allocate(Uwrk0(1:rk), source=Uwrk0(1)); call mat_zero(Uwrk0)
         if (.not. allocated(Swrk0)) allocate(Swrk0(1:rk,1:rk)); Swrk0 = 0.0_wp
         call M_forward_map        (U,     S,     LTI,          0.5*tau, info, exptA, iftrans=.true.)
         ! Predictor step
         call G_forward_map_riccati(Uwrk0, Swrk0, LTI, Q, Rinv,     tau, info)
         ! Second order integration
         call G_forward_map_riccati(U,     S,     LTI, Q, Rinv,     tau, info, Uwrk0, Swrk0)
         call M_forward_map        (U,     S,     LTI,          0.5*tau, info, exptA, iftrans=.true.)
      end select

      return

   end subroutine numerical_low_rank_splitting_riccati_step

   subroutine G_forward_map_riccati(U, S, LTI, Q, Rinv, tau, info, Upred, Spred)
      !> Low-rank factors
      class(abstract_vector),           intent(inout) :: U(:)   ! basis
      real(kind=wp),                    intent(inout) :: S(:,:) ! coefficients
      ! LTI system
      class(abstract_lti_system),       intent(in)    :: LTI
      class(abstract_linop),            intent(in)    :: Q
      real(kind=wp),                    intent(in)    :: Rinv(:,:)
      !> Integration step size         
      real(kind=wp),                    intent(in)    :: tau
      !> Information flag
      integer,                          intent(out)   :: info
      !> Intermediate values
      class(abstract_vector), optional, intent(in)    :: Upred(:)   ! basis
      class(abstract_vector),           allocatable   :: Up(:)  
      real(kind=wp),          optional, intent(in)    :: Spred(:,:) ! coefficients
      real(kind=wp),                    allocatable   :: Sp(:,:)
      
      !> Local variables
      class(abstract_vector),           allocatable   :: U1(:)
      class(abstract_vector),           allocatable   :: QU(:)
      real(kind=wp),                    allocatable   :: UTB(:,:)
      integer                                         :: m, rk

      rk = size(U)
      m  = size(LTI%B)
      if (.not. allocated(U1))  allocate(U1(1:rk), source=U(1)); call mat_zero(U1)
      if (.not. allocated(QU))  allocate(QU(1:rk), source=U(1)); call mat_zero(QU) 
      if (.not. allocated(UTB)) allocate(UTB(1:rk,1:m)); UTB = 0.0_wp

      if ( present(Upred) .and. present(Spred) ) then
         ! second order in time
         allocate(Up(1:rk), source=Upred(1:rk))
         allocate(Sp(1:rk,1:rk)); Sp = Spred
         ! Compute the intermediate state at t = t0 + dt/2 as the 
         ! average of states at t0 (U,S) and t0 + dt (Up,Sp)
         call mat_axpby(Up, 0.5_wp, U, 0.5_wp)
         call mat_axpby(Sp, 0.5_wp, S, 0.5_wp)
         ! first steps based on previous state (U,S)
         call K_step_riccati(U1, S,  QU, UTB, U,      LTI, Q, Rinv, 0.5*tau, info)
         call S_step_riccati(    S,  QU, UTB, U,  U1, LTI, Q, Rinv, 0.5*tau, info)
         call L_step_riccati(    S,           U,  U1, LTI, Q, Rinv,     tau, info)
         ! second set in reverse order based on intermediate state (Up,Sp)
         call S_step_riccati(    Sp, QU, UTB, Up, U1, LTI, Q, Rinv, 0.5*tau, info, reverse=.true.)
         call K_step_riccati(U1, Sp, QU, UTB, Up,     LTI, Q, Rinv, 0.5*tau, info, reverse=.true.)
      else
         ! first order in time
         call K_step_riccati(U1, S,  QU, UTB, U,      LTI, Q, Rinv,     tau, info)
         call S_step_riccati(    S,  QU, UTB, U,  U1, LTI, Q, Rinv,     tau, info)
         call L_step_riccati(    S,           U,  U1, LTI, Q, Rinv,     tau, info)
      end if
      
      !> Copy data to output
      call mat_copy(U, U1)
               
      return
   end subroutine G_forward_map_riccati

   subroutine K_step_riccati(U1, S, QU, UTB, U, LTI, Q, Rinv, tau, info, reverse)
      !> Low-rank factors
      class(abstract_vector),     intent(out)   :: U1(:)  ! basis
      real(kind=wp),              intent(inout) :: S(:,:) ! coefficients
      class(abstract_vector),     intent(inout) :: QU(:)  ! basis
      real(kind=wp),              intent(inout) :: UTB(:,:) ! matrix
      !> Low-rank factors
      class(abstract_vector),     intent(in)    :: U(:)   ! basis
      ! LTI system
      class(abstract_lti_system), intent(in)    :: LTI
      class(abstract_linop),      intent(in)    :: Q
      real(kind=wp),              intent(in)    :: Rinv(:,:)
      !> Integration step size
      real(kind=wp),              intent(in)    :: tau
      !> Information flag
      integer,                    intent(out)   :: info
      !> in second order scheme, are we in the forward or reverse branch?
      logical, optional,          intent(in)    :: reverse
      logical                                   :: reverse_order

      !> Local variables
      class(abstract_vector),     allocatable   :: Uwrk0(:)
      real(kind=wp),              allocatable   :: Swrk0(:,:)
      real(kind=wp),              allocatable   :: BTU(:,:)
      real(kind=wp),              allocatable   :: PTwrk(:,:)
      real(kind=wp),              allocatable   :: Qwrk(:,:)
      integer,                    allocatable   :: perm(:)   ! Permutation vector
      integer                                   :: i, m, rk

      ! Optional arguments
      reverse_order = optval(reverse, .false.)

      info = 0

      rk = size(U)
      m  = size(LTI%B)
      if (.not. allocated(Uwrk0)) allocate(Uwrk0(1:rk), source=U(1)); call mat_zero(Uwrk0)
      if (.not. allocated(Swrk0)) allocate(Swrk0(1:rk,1:rk)); Swrk0 = 0.0_wp
      if (.not. allocated(BTU))   allocate(BTU(  1:m, 1:rk)); BTU   = 0.0_wp
      if (.not. allocated(perm))  allocate(perm( 1:rk));      perm  = 0

      call mat_mult(U1, U, S)                         ! K0

      ! Constant part --> QU
      if (.not. reverse_order) then
         ! compute QU and UTB and pass to S step
         call mat_zero(QU)
         do i = 1, rk
            call Q%matvec(U(i), QU(i))
         end do
         call mat_mult(UTB, U, LTI%B)                                      !              U0.T @ B
      end if

      ! non-linear part --> Uwrk0
      call mat_mult(BTU, LTI%B, U1)                                        !                                      B.T @ U0
      call mat_mult(Uwrk0, U1, matmul(UTB, matmul(Rinv, matmul(BTU, S))))  ! (U0 @ S0) @ (U0.T @ B) @ (R^(-1) @ ((B.T @ U0) @ S0))
      
      ! combine to form G( K @ U.T ) @ U --> Uwrk0
      call mat_axpby(Uwrk0, -1.0_wp, QU, 1.0_wp)

      !> Construct solution U1
      call mat_axpby(U1, 1.0_wp, Uwrk0, tau) ! K0 + tau*Kdot

      !> Orthonormalize in-place
      call qr_factorization(U1, Swrk0, perm, info, ifpivot = .true.)
      call apply_permutation(Swrk0, perm, trans = .true.)
      S = Swrk0

      return
   end subroutine K_step_riccati

   subroutine S_step_riccati(S, QU, UTB, U, U1, LTI, Q, Rinv, tau, info, reverse)
      !> Low-rank factors
      real(kind=wp),              intent(inout) :: S(:,:) ! coefficients
      class(abstract_vector),     intent(inout) :: QU(:)  ! basis
      real(kind=wp),              intent(inout) :: UTB(:,:) ! matrix
      !> Low-rank factors
      class(abstract_vector),     intent(in)    :: U(:)   ! old basis
      class(abstract_vector),     intent(in)    :: U1(:)  ! updated basis
      ! LTI system
      class(abstract_lti_system), intent(in)    :: LTI
      class(abstract_linop),      intent(in)    :: Q
      real(kind=wp),              intent(in)    :: Rinv(:,:)
      !> Integration step size
      real(kind=wp),              intent(in)    :: tau
      !> Information flag
      integer,                    intent(out)   :: info
      !> in second order scheme, are we in the forward or reverse branch?
      logical, optional,          intent(in)    :: reverse
      logical                                   :: reverse_order

      !> Local variables
      real(kind=wp),              allocatable   :: BTU(:,:)
      real(kind=wp),              allocatable   :: Pwrk1(:,:)
      real(kind=wp),              allocatable   :: PTwrk(:,:)
      real(kind=wp),              allocatable   :: Swrk0(:,:)
      real(kind=wp),              allocatable   :: Swrk1(:,:)
      integer                                   :: i, m, rk

      info = 0

      ! Optional arguments
      reverse_order = optval(reverse, .false.)

      rk = size(U)
      m  = size(LTI%B)
      if (.not. allocated(BTU))   allocate(BTU  (1:m, 1:rk)); BTU   = 0.0_wp
      if (.not. allocated(Swrk0)) allocate(Swrk0(1:rk,1:rk)); Swrk0 = 0.0_wp
      if (.not. allocated(Swrk1)) allocate(Swrk1(1:rk,1:rk)); Swrk1 = 0.0_wp 

      ! Constant part --> Swrk0
      if (reverse_order) then
         ! compute QU and UTB pass it to K step
         call mat_zero(QU)
         do i = 1, rk
            call Q%matvec(U(i), QU(i))
         end do
         call mat_mult(UTB, U, LTI%B)                              !       U0.T @ B
      endif
      call mat_mult(Swrk0, U1, QU)

      ! non-linear part --> Swrk1
      call mat_mult(BTU, LTI%B, U1)                                !                              B.T @ U1
      Swrk1 = matmul(S, matmul(UTB, matmul(Rinv, matmul(BTU, S)))) ! S0 @ (U0.T @ B @ (R^(-1) @ ((B.T @ U1) @ S0)))

      ! combine to form -U1.T @ G( U1 @ S @ U0.T ) @ U0
      call mat_axpby(Swrk0, -1.0_wp, Swrk1, 1.0_wp)

      !> Construct solution S1
      call mat_axpby(S, 1.0_wp, Swrk0, tau)

      return
   end subroutine S_step_riccati

   subroutine L_step_riccati(S, U, U1, LTI, Q, Rinv, tau, info)
      !> Low-rank factors
      real(kind=wp),              intent(inout) :: S(:,:) ! coefficients
      !> Low-rank factors
      class(abstract_vector),     intent(inout) :: U(:)   ! basis before Kstep
      class(abstract_vector),     intent(in)    :: U1(:)   ! basis after Kstep
      ! LTI system
      class(abstract_lti_system), intent(in)    :: LTI
      class(abstract_linop),      intent(in)    :: Q
      real(kind=wp),              intent(in)    :: Rinv(:,:)
      !> Integration step size
      real(kind=wp),              intent(in)    :: tau
      !> Information flag
      integer,                    intent(out)   :: info

      !> Local variables
      class(abstract_vector),     allocatable   :: Uwrk0(:)
      class(abstract_vector),     allocatable   :: Uwrk1(:)
      class(abstract_vector),     allocatable   :: Uwrk2(:)
      real(kind=wp),              allocatable   :: Swrk0(:,:)
      real(kind=wp),              allocatable   :: Swrk1(:,:)
      real(kind=wp),              allocatable   :: Qwrk(:,:)
      integer                                   :: i, m, rk

      info = 0

      rk = size(U)
      m  = size(LTI%B)
      if (.not. allocated(Uwrk0)) allocate(Uwrk0(1:rk), source=U(1)); call mat_zero(Uwrk0)
      if (.not. allocated(Uwrk1)) allocate(Uwrk1(1:rk), source=U(1)); call mat_zero(Uwrk1)
      if (.not. allocated(BTU))   allocate(BTU(  1:m, 1:rk)); BTU   = 0.0_wp
      if (.not. allocated(UTB))   allocate(UTB(  1:rk,1:m )); UTB   = 0.0_wp
      if (.not. allocated(Swrk0)) allocate(Swrk0(1:rk,1:rk)); Swrk0 = 0.0_wp

      call mat_mult(Uwrk1, U, transpose(S)) ! L0.T = U0 @ S.T

      ! Constant part --> Uwrk0
      call mat_zero(Uwrk0)
      do i = 1, rk
         call Q%matvec(U1(i), Uwrk0(i))
      end do

      ! non-linear part --> U
      call mat_mult(BTU, LTI%B, Uwrk1)       !                                     B.T @ U0 @ S.T
      call mat_mult(UTB, U1, LTI%B)          !               U1.T @ B
      Swrk0 = matmul(UTB, matmul(Rinv, BTU)) !              (U1.T @ B) @ (R^(-1) @ (B.T @ U0 @ S.T))
      call mat_mult(U, Uwrk1, Swrk0)         ! (U0 @ S.T) @ (U1.T @ B  @  R^(-1) @  B.T @ U0 @ S.T)

      ! combine to form U1.T @ G( U1.T@L.T )
      call mat_axpby(Uwrk0, 1.0_wp, U, -1.0_wp)

      !> Construct solution U1
      call mat_axpby(Uwrk1, 1.0_wp, Uwrk0, tau)       ! L0 + tau*Ldot

      !> Update S
      call mat_mult(S, Uwrk1, U1)

      return
   end subroutine L_step_riccati

end module lightROM_RiccatiSolvers
