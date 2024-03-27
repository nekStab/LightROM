module LightROM_RiccatiSolvers
   use LightKrylov
   use LightKrylov_expmlib
   use lightkrylov_BaseKrylov
   use LightROM_LyapunovUtils
   use LightROM_LyapunovSolvers
   use LightROM_utils
   use stdlib_optval, only : optval
   use stdlib_linalg, only : eye
   implicit none

   !> work arrays
   class(abstract_vector),  allocatable   :: Uwrk0(:)
   class(abstract_vector),  allocatable   :: Uwrk1(:)
   class(abstract_vector) , allocatable   :: Uwrk2(:)
   class(abstract_vector),  allocatable   :: CTQCU(:)
   real(kind=wp),           allocatable   :: Swrk0(:,:)
   real(kind=wp),           allocatable   :: Swrk1(:,:)
   real(kind=wp),           allocatable   :: Swrk2(:,:)

   private
   public :: numerical_low_rank_splitting_integrator_riccati

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
   ! - Ctrans : Low-rank input for the Riccati equation   [Input]
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
   subroutine numerical_low_rank_splitting_integrator_riccati(U,S,A,B,Ctrans,Tend,tau,torder,info,exptA)
      !> Low-rank factors
      class(abstract_vector),              intent(inout) :: U(:) ! basis
      real(kind=wp),                       intent(inout) :: S(:,:) ! coefficients
      !> System
      class(abstract_linop),               intent(in)    :: A    ! linear operator
      class(abstract_vector),              intent(in)    :: B(:) !
      class(abstract_vector) ,             intent(in)    :: Ctrans(:)
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
      integer                                :: rk, istep, nsteps, iostep
      real(kind=wp)                          :: T
      logical, parameter                     :: verbose = .false.
      procedure(abstract_exptA), pointer     :: p_exptA => null()

      T = 0.0_wp
      rk = size(U)

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
         call numerical_low_rank_splitting_step(U, S, A, B, Ctrans, tau, torder, info, p_exptA)

         T = T + tau
         !> here we should do some checks such as whether we have reached steady state
         if ( mod(istep,iostep) .eq. 0 ) then
            if (verbose) then
               write(*, *) "INFO : ", ISTEP, " steps of DLRA computed. T = ",T
            endif
         endif
      enddo dlra
      return

   end subroutine numerical_low_rank_splitting_integrator_riccati

   !-----------------------------
   !-----     UTILITIES     -----
   !-----------------------------
   subroutine numerical_low_rank_splitting_riccati_step(U, S, A, B, Ctrans, tau, torder, info, exptA)
      !> Low-rank factors
      class(abstract_vector),    intent(inout) :: U(:) ! basis
      real(kind=wp),             intent(inout) :: S(:,:) ! coefficients
      !> Linear operator and rhs Cholesky factor
      class(abstract_linop),     intent(in)    :: A    ! linear operator
      class(abstract_vector),    intent(in)    :: B(:)
      class(abstract_vector) ,   intent(in)    :: Ctrans(:)
      !> Integration step size
      real(kind=wp),             intent(in)    :: tau
      !> Order of time integration
      integer,                   intent(in)    :: torder
      !> Information flag
      integer,                   intent(out)   :: info
      !> Routine for computation of the exponential propagator
      procedure(abstract_exptA)                :: exptA

      !> Local variables
      integer        ::   istep, nsteps, integrator
      class(abstract_vector),    allocatable   :: Uwrk0(:) ! basis
      real(kind=wp),             allocatable   :: Swrk0(:,:) ! coefficients

      if (.not. allocated(Uwrk0)) allocate(Uwrk0(1:rk), source=Uwrk0(1)); 
      if (.not. allocated(Swrk0)) allocate(Swrk0(1:rk,1:rk));
      call mat_zero(Uwrk0); Swrk0 = 0.0_wp

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
         call M_forward_map_riccati(U, S, A, tau, info, exptA)
         call G_forward_map_riccati(U, S, B, Ctrans, tau, info)
      case (2) ! Strang splitting
         call M_forward_map_riccati(U, S, A, 0.5*tau, info, exptA)
         ! Predictor step
         call G_forward_map_riccati(Uwrk0, Swrk0, B, Ctrans, tau, info)
         ! Second order integration
         call G_forward_map_riccati(U, S, B, Ctrans, tau, Up=Uwrk0, Sp=Swrk0, info)
         call M_forward_map_riccati(U, S, A, 0.5*tau, info, exptA)
      end select

      return

   end subroutine numerical_low_rank_splitting_riccati_step

   subroutine G_forward_map_riccati(U, S, B, Ctrans, tau, Up, Sp, info)
      !> Low-rank factors
      class(abstract_vector),           intent(inout) :: U(:)   ! basis
      real(kind=wp),                    intent(inout) :: S(:,:) ! coefficients
      !> low-rank inputs         
      class(abstract_vector),           intent(in)    :: B(:)
      class(abstract_vector),           intent(in)    :: Ctrans(:)
      !> Integration step size         
      real(kind=wp),                    intent(in)    :: tau
      class(abstract_vector), optional, intent(in)    :: Up(:)   ! basis
      real(kind=wp),          optional, intent(in)    :: Sp(:,:) ! coefficients
      !> Information flag
      integer,                          intent(out)   :: info
   
      !> Local variables
      class(abstract_vector),           allocatable   :: U1(:)
      class(abstract_vector),           allocatable   :: CTQCU(:)
      integer                                         :: rk

         torder = 2
      else
         torder = 1
      end if   

      rk = size(U)
      if (.not. allocated(U1))    allocate(U1(1:rk),    source=U(1)); 
      if (.not. allocated(CTQCU)) allocate(CTQCU(1:rk), source=U(1));
      call mat_zero(U1); call mat_zero(CTQCU)

      if ( present(Up) .and. present(Sp) ) then
         ! second order in time
         ! Compute the intermediate state at t = t0 + dt/2 as the 
         ! average of states at t0 (U,S) and t0 + dt (Up,Sp)
         call mat_axpby(Up, 0.5_wp, U, 0.5_wp)
         call mat_axpby(Sp, 0.5_wp, S, 0.5_wp)
         ! first steps based on previous state (U,S)
         call K_step(U1, S,  CTQCU, U,      B, 0.5*tau, info)
         call S_step(    S,  CTQCU, U,  U1, B, 0.5*tau, info)
         call L_step(    S,         U,  U1, B,     tau, info)
         ! second set in reverse order based on intermediate state (Up,Sp)
         call S_step(    Sp, CTQCU, Up, U1, B, 0.5*tau, reverse = .true., info)
         call K_step(U1, Sp, CTQCU, Up,     B, 0.5*tau, reverse = .true., info)
      else
         ! first order in time
         call K_step(U1, S, CTQCU, U,     B, tau, info)
         call S_step(    S, CTQCU, U, U1, B, tau, info)
         call L_step(    S,        U, U1, B, tau, info)
      end if
      
      !> Copy data to output
      call mat_copy(U, U1)
               
      return
   end subroutine G_forward_map_riccati

   subroutine K_step_riccati(U1, S, CTQCU, U, B, Ctrans, tau, reverse, info)
      !> Low-rank factors
      class(abstract_vector) , intent(out)   :: U1(:)  ! basis
      real(kind=wp)          , intent(inout) :: S(:,:) ! coefficients
      class(abstract_vector) , intent(inout) :: CTQCU(:)  ! basis
      !> Low-rank factors
      class(abstract_vector) , intent(in)    :: U(:)   ! basis
      !> low-rank inputs
      class(abstract_vector) , intent(in)    :: B(:)
      class(abstract_vector) , intent(in)    :: Ctrans(:)
      !> Integration step size
      real(kind=wp)          , intent(in)    :: tau
      !> Information flag
      integer                , intent(out)   :: info
      !> in second order scheme, are we in the forward or reverse branch?
      logical, optional      , intent(in)    :: reverse
      logical,                               :: reverse_order

      !> Local variables
      class(abstract_vector) , allocatable   :: Uwrk0(:)
      class(abstract_vector) , allocatable   :: Uwrk1(:)
      real(kind=wp)          , allocatable   :: Swrk0(:,:)
      real(kind=wp)          , allocatable   :: Swrk1(:,:)
      integer                , allocatable   :: perm(:)   ! Permutation vector
      integer                                :: rk

      ! Optional arguments
      reverse_order = optval(reverse, .false.)

      info = 0

      rk = size(U)
      if (.not. allocated(Uwrk0)) allocate(Uwrk0(1:rk), source=U(1));
      if (.not. allocated(Uwrk1)) allocate(Uwrk1(1:rk), source=U(1));
      if (.not. allocated(Swrk0)) allocate(Swrk0(1:rk,1:rk)); 
      if (.not. allocated(Swrk1)) allocate(Swrk1(1:rk,1:rk)); 
      allocate(perm(1:rk)); perm = 0
      call mat_zero(Uwrk0); call mat_zero(Uwrk1)
      Swrk0 = 0.0_wp; Swrk1 = 0.0_wp

      call mat_mult(U1, U, S)                         ! K0

      ! Constant part --> Uwrk0
      if (.not. reverse_order) then
         ! compute CTQCU and pass it to S step
         call mat_mult(Uwrk1, Ctrans, U)                 ! C @ U0
         call apply_Q(Uwrk1)                             ! Q @ C @ U0             (in-place, user-defined)
         ! apply_Q could also be a function (in-place, user-defined)
         call mat_mult(CTQCU, Ctrans, Uwrk1)             ! C.T @ Q @ C @ U0
      end if

      ! non-linear part --> Uwrk1
      call mat_mult(Swrk0, B, U1)                     ! B.T @ U0 @ S0
      call apply_Rinv(Swrk0)                          ! R^(-1) @ B.T @ U0 @ S0    (in-place, user-defined)
      ! apply_Rinv could also be a function (in-place, user-defined)
      call mat_mult(Swrk1, U, B)                      ! U0.T @ B
      call mat_mult(Uwrk1, U1, matmul(Swrk0, Swrk1))  ! (U0 @ S0) @ (U0.T @ B) @ (R^(-1) @ B.T @ U0 @ S0)
      
      ! combine to form G( K @ U.T ) @U --> Uwrk1
      call mat_axpby(Uwrk1, -1.0_wp, CTQCU, 1.0_wp)

      !> Construct solution U1
      call mat_axpby(U1, 1.0_wp, Uwrk1, tau) ! K0 + tau*Kdot

      !> Orthonormalize in-place
      call qr_factorization(U1, Swrk0, perm, info, ifpivot = .true.)
      call apply_permutation(Swrk0, perm, trans = .true.)
      S = Swrk0

      return
   end subroutine K_step_riccati

   subroutine S_step_riccati(S, CTQCU, U, U1, B, tau, reverse, info)
      !> Low-rank factors
      real(kind=wp)          , intent(inout) :: S(:,:) ! coefficients
      class(abstract_vector) , intent(inout) :: CTQCU(:)  ! basis
      !> Low-rank factors
      class(abstract_vector) , intent(in)    :: U(:)   ! old basis
      class(abstract_vector) , intent(in)    :: U1(:)  ! updated basis
      !> low-rank inputs
      class(abstract_vector) , intent(in)    :: B(:)
      !> Integration step size
      real(kind=wp)          , intent(in)    :: tau
      !> Information flag
      integer                , intent(out)   :: info
      !> in second order scheme, are we in the forward or reverse branch?
      logical, optional      , intent(in)    :: reverse
      logical,                               :: reverse_order

      !> Local variables
      real(kind=wp)          , allocatable   :: Swrk0(:,:)
      real(kind=wp)          , allocatable   :: Swrk1(:,:)
      real(kind=wp)          , allocatable   :: Swrk2(:,:)
      integer                                :: rk

      info = 0

      ! Optional arguments
      reverse_order = optval(reverse, .false.)

      rk = size(U)
      if (.not. allocated(Swrk0)) allocate(Swrk0(1:rk,1:rk)); 
      if (.not. allocated(Swrk1)) allocate(Swrk1(1:rk,1:rk)); 
      if (.not. allocated(Swrk2)) allocate(Swrk2(1:rk,1:rk)); 
      Swrk0 = 0.0_wp; Swrk1 = 0.0_wp; Swrk2 = 0.0_wp

      ! non-linear part --> Swrk1
      call mat_mult(Swrk0, U, B)                      ! U0.T @ B
      call mat_mult(Swrk1, B, U1)                     ! B.T @ U1
      Swrk2 = matmul(Swrk1, S)
      call apply_Rinv(Swrk2)                          ! R^(-1) @ B.T @ U1 @ S0    (in-place, user-defined)
      ! apply_Rinv could also be a function (in-place, user-defined)
      Swrk1 = matmul(S, matmul(Swrk0, Swrk2))         ! S0 @ ( U0.T @ B ) @ ( R^(-1) @ B.T @ U1 @ S0 ) 

      ! Constant part --> Swrk0
      if (reverse_order) then
         ! compute CTQCU and pass it to K step
         call mat_mult(Uwrk1, Ctrans, U)                 ! C @ U0
         call apply_Q(Uwrk1)                             ! Q @ C @ U0             (in-place, user-defined)
         ! apply_Q could also be a function (in-place, user-defined)
         call mat_mult(CTQCU, Ctrans, Uwrk1)             ! C.T @ Q @ C @ U0
      endif
      call mat_mult(Swrk0, U1, CTQCU)                    ! U1.T @ C.T @ Q @ C @ U0

      ! combine to form -U1.T @ G( U1 @ S @ U0.T ) @ U0
      call mat_axpby(Swrk0, -1.0_wp, Swrk1, 1.0_wp)

      !> Construct solution S1
      call mat_axpby(S, 1.0_wp, Swrk, tau)

      return
   end subroutine S_step_riccati

   subroutine L_step_riccati(S, U, U1, B, Ctrans, tau, info)
      !> Low-rank factors
      real(kind=wp)          , intent(inout) :: S(:,:) ! coefficients
      !> Low-rank factors
      class(abstract_vector) , intent(inout) :: U(:)   ! basis before Kstep
      class(abstract_vector) , intent(in)    :: U1(:)   ! basis after Kstep
      !> low-rank inputs
      class(abstract_vector) , intent(in)    :: B(:)
      class(abstract_vector) , intent(in)    :: Ctrans(:)
      !> Integration step size
      real(kind=wp)          , intent(in)    :: tau
      !> Information flag
      integer                , intent(out)   :: info

      !> Local variables
      class(abstract_vector) , allocatable   :: Uwrk0(:)
      class(abstract_vector) , allocatable   :: Uwrk1(:)
      class(abstract_vector) , allocatable   :: Uwrk2(:)
      real(kind=wp)          , allocatable   :: Swrk0(:,:)
      real(kind=wp)          , allocatable   :: Swrk1(:,:)
      integer                                :: rk

      info = 0

      rk = size(U)
      if (.not. allocated(Uwrk0)) allocate(Uwrk0(1:rk), source=U(1));
      if (.not. allocated(Uwrk1)) allocate(Uwrk1(1:rk), source=U(1));
      if (.not. allocated(Uwrk2)) allocate(Uwrk2(1:rk), source=U(1));
      if (.not. allocated(Swrk0)) allocate(Swrk0(1:rk,1:rk)); 
      if (.not. allocated(Swrk1)) allocate(Swrk1(1:rk,1:rk)); 
      call mat_zero(Uwrk0); call mat_zero(Uwrk1); call mat_zero(Uwrk2)
      Swrk0 = 0.0_wp; Swrk1 = 0.0_wp;

      ! Constant part --> Uwrk0
      call mat_mult(Uwrk1, Ctrans, U1)                ! C @ U1
      call apply_Q(Uwrk1)                             ! Q @ C @ U1                 (in-place, user-defined)
      ! apply_Q could also be a function (in-place, user-defined)
      call mat_mult(Uwrk0, Ctrans, Uwrk1)             ! C.T @ Q @ C @ U1

      ! non-linear part --> Uwrk1
      call mat_mult(Uwrk1, U, transpose(S))           ! L0.T
      call mat_mult(Swrk0, B, Uwrk1)                  ! B.T @ U0 @ S.T
      call apply_Rinv(Swrk0)                          ! R^(-1) @ B.T @ U0 @ S.T    (in-place, user-defined)
      ! apply_Rinv could also be a function (in-place, user-defined)
      call mat_mult(Swrk1, U1, B)                     ! U1.T @ B
      call mat_mult(U, Uwrk1, Swrk1)                  ! (U0 @ S.T) @ (U1.T @ B)
      call mat_mult(Uwrk2, U, Swrk0)                  ! (U0 @ S.T) @ (U1.T @ B) @ (R^(-1) @ B.T @ U0 @ S.T)

      ! combine to form U1.T @ G( U1.T@L.T )
      call mat_axpby(Uwrk0, 1.0_wp, Uwrk2, -1.0_wp)

      !> Construct solution U1
      call mat_axpby(Uwrk1, 1.0_wp, Uwrk0, tau)       ! L0 + tau*Ldot

      !> Update S
      call mat_mult(S, Uwrk1, U1)

      return
   end subroutine L_step_riccati

end module lightROM_RiccatiSolvers
