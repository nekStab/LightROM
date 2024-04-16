module LightROM_RiccatiSolvers
   use LightKrylov
   use LightKrylov_utils, only : print_mat
   use LightKrylov_expmlib
   use lightkrylov_BaseKrylov
   use LightROM_LyapunovUtils
   use LightROM_LyapunovSolvers
   use LightROM_RiccatiUtils
   !use LightROM_Utils
   use lightrom_AbstractLTIsystems
   use stdlib_optval, only : optval
   use stdlib_linalg, only : eye
   implicit none

   !> work arrays
   class(abstract_vector),  allocatable   :: Uwrk0(:)
   class(abstract_vector),  allocatable   :: Uwrk1(:)
   class(abstract_vector),  allocatable   :: U1(:)
   class(abstract_vector),  allocatable   :: QU(:)
   real(kind=wp),           allocatable   :: Swrk0(:,:)
   real(kind=wp),           allocatable   :: Swrk1(:,:)
   ! prediction step
   class(abstract_vector),  allocatable   :: U0(:)
   class(abstract_vector),  allocatable   :: T0(:)
   class(abstract_vector),  allocatable   :: Ut(:)
   class(abstract_vector),  allocatable   :: Tt(:)
   real(kind=wp),           allocatable   :: S0(:,:)

   private
   public :: numerical_low_rank_splitting_riccati_integrator
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
   subroutine numerical_low_rank_splitting_riccati_integrator(X,LTI,Qc,Rinv,Tend,tau,torder,info,exptA,iftrans)
      !> Low-rank state
      class(abstract_sym_low_rank_state),  intent(inout) :: X
      ! LTI system
      class(abstract_lti_system),          intent(in)    :: LTI
      real(kind=wp),                       intent(in)    :: Qc(:,:)
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
      !> use transpose
      logical,                   optional, intent(in)    :: iftrans
      logical                                            :: trans

      !> Local variables   
      integer                                            :: istep, nsteps, iostep
      real(kind=wp)                                      :: T
      logical, parameter                                 :: verbose = .false.
      procedure(abstract_exptA), pointer                 :: p_exptA => null()

      ! Optional argument
      trans = optval(iftrans, .false.)

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
         call numerical_low_rank_splitting_riccati_step(X, LTI, Qc, Rinv, tau, torder, info, p_exptA, trans)

         T = T + tau
         !> here we should do some checks such as whether we have reached steady state
         if ( mod(istep,iostep) .eq. 0 ) then
            if (verbose) then
               write(*, *) "INFO : ", ISTEP, " steps of DLRA computed. T = ",T
            endif
         endif
      enddo dlra

      if (allocated(Uwrk0)) deallocate(Uwrk0)
      if (allocated(Uwrk1)) deallocate(Uwrk1)
      if (allocated(Swrk0)) deallocate(Swrk0)
      if (allocated(Swrk1)) deallocate(Swrk1)
      if (allocated(QU)) deallocate(QU)
      if (allocated(U1)) deallocate(U1)
      if (allocated(U0)) deallocate(U0)
      if (allocated(T0)) deallocate(T0)
      if (allocated(Ut)) deallocate(Ut)
      if (allocated(Tt)) deallocate(Tt)
      if (allocated(S0)) deallocate(S0)

      return

   end subroutine numerical_low_rank_splitting_riccati_integrator

   !-----------------------------
   !-----     UTILITIES     -----
   !-----------------------------
   subroutine numerical_low_rank_splitting_riccati_step(X, LTI, Qc, Rinv, tau, torder, info, exptA, iftrans)
      !> Low-rank state
      class(abstract_sym_low_rank_state),  intent(inout) :: X
      ! LTI system
      class(abstract_lti_system),          intent(in)    :: LTI
      real(kind=wp),                       intent(in)    :: Qc(:,:)
      real(kind=wp),                       intent(in)    :: Rinv(:,:)
      !> Integration step size
      real(kind=wp),                       intent(in)    :: tau
      !> Order of time integration
      integer,                             intent(in)    :: torder
      !> Information flag
      integer,                             intent(out)   :: info
      !> Routine for computation of the exponential propagator
      procedure(abstract_exptA)                          :: exptA
      !> use transpose
      logical,                   optional, intent(in)    :: iftrans
      logical                                            :: trans

      !> Local variables
      integer                                            :: istep, nsteps, integrator, rk

      ! Optional argument
      trans = optval(iftrans, .false.)

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
         call M_forward_map        (X, LTI,               tau, info, exptA, trans)
         call G_forward_map_riccati(X, LTI, Qc, Rinv,     tau, info)
      case (2) ! Strang splitting
         ! Prepare arrays for predictor step
         rk = size(X%U)
         ! precomputation arrays
         if (.not. allocated(U0)) allocate(U0(1:rk), source=X%U(1:rk))
         if (.not. allocated(T0)) allocate(T0(1:rk), source=X%U(1:rk))
         if (.not. allocated(Ut)) allocate(Ut(1:rk), source=X%U(1:rk))
         if (.not. allocated(Tt)) allocate(Tt(1:rk), source=X%U(1:rk))
         if (.not. allocated(S0)) allocate(S0(1:rk,1:rk))
         call mat_zero(U0); call mat_zero(T0); call mat_zero(Ut); call mat_zero(Tt); S0 = 0.0_wp
         ! scratch arrays
         if (.not. allocated(Uwrk0)) allocate(Uwrk0(1:rk), source=X%U(1:rk))
         if (.not. allocated(Swrk0)) allocate(Swrk0(1:rk,1:rk))
         call mat_zero(Uwrk0); Swrk0 = 0.0_wp
         ! second order step
         call M_forward_map        (X, LTI,           0.5*tau, info, exptA, trans)
         ! Save current state
         call mat_copy(U0, X%U); S0 = X%S                             ! --> save      
         ! Precompute T0
         call mat_mult(Uwrk0, U0, S0)                                 ! K0 = U0 @ S0
         call apply_p_outerproduct_w(Swrk0, U0, Uwrk0, LTI%B, Rinv)   ! (U0.T) @ B @ R^(-1) @ B.T @ K0
         call mat_mult(T0, Uwrk0, Swrk0)                              ! K0 @ Swrk0
         ! First order integration
         call G_forward_map_riccati(X, LTI, Qc, Rinv,     tau, info, ifpred=.true., T0=T0)
         ! Precompute Tt
         call mat_copy(Ut, X%U)                                       ! --> save
         call mat_mult(Uwrk0, X%U, X%S)                               ! Kt = Ut @ St
         call apply_p_outerproduct_w(Swrk0, X%U, Uwrk0, LTI%B, Rinv)  ! (Ut.T) @ B @ R^(-1) @ B.T @ Kt
         call mat_mult(Tt, Uwrk0, Swrk0)                              ! Kt @ Swrk0
         ! Reset state
         call mat_copy(X%U, U0); X%S = S0
         ! Second order integration
         call G_forward_map_riccati(X, LTI, Qc, Rinv,     tau, info, ifpred=.false., T0=T0, Tt=Tt, U0=U0, Ut=Ut)
   !      call G_forward_map_riccati(X, LTI, Qc, Rinv,     tau, info)
         call M_forward_map        (X, LTI,           0.5*tau, info, exptA, trans)
      end select

      return

   end subroutine numerical_low_rank_splitting_riccati_step

   subroutine G_forward_map_riccati(X, LTI, Qc, Rinv, tau, info, ifpred, T0, Tt, U0, Ut)
      !> Low-rank state
      class(abstract_sym_low_rank_state),  intent(inout) :: X      ! current state
      ! LTI system
      class(abstract_lti_system),          intent(in)    :: LTI
      real(kind=wp),                       intent(in)    :: Qc(:,:)
      real(kind=wp),                       intent(in)    :: Rinv(:,:)
      !> Integration step size         
      real(kind=wp),                       intent(in)    :: tau
      !> Information flag
      integer,                             intent(out)   :: info
      !> for second order: predictor or corrector step?
      logical,                optional,    intent(in)    :: ifpred
      !> Intermediate values
      class(abstract_vector), optional,    intent(inout) :: T0(:)  ! will be reused as Gamma
      class(abstract_vector), optional,    intent(in)    :: Tt(:)
      class(abstract_vector), optional,    intent(in)    :: U0(:)
      class(abstract_vector), optional,    intent(in)    :: Ut(:)
      
      !> Local variables
      integer                                            :: rk

      rk = size(X%U)
      if (.not. allocated(U1))  allocate(U1( 1:rk), source=X%U(1)); 
      if (.not. allocated(QU))  allocate(QU( 1:rk), source=X%U(1));
      if (.not. allocated(Swrk0)) allocate(Swrk0(1:rk,1:rk))
      call mat_zero(U1); call mat_zero(QU)

      if (present(ifpred)) then
         ! second order in time
         if (ifpred) then ! predictor step with precomputed T0
            call K_step_riccati(X, U1, QU, LTI, Qc, Rinv,     tau, info, reverse=.false., NL=T0)
            call S_step_riccati(X, U1, QU, LTI, Qc, Rinv,     tau, info, reverse=.false.)
            call L_step_riccati(X, U1,     LTI, Qc, Rinv,     tau, info)
         else             ! corrector step with precomputed T0, Tt and U0, Ut
            ! forward steps based on T0, U0 (within X)
            call K_step_riccati(X, U1, QU, LTI, Qc, Rinv, 0.5*tau, info, reverse=.false., NL=T0)
            call S_step_riccati(X, U1, QU, LTI, Qc, Rinv, 0.5*tau, info, reverse=.false.)
            call L_step_riccati(X, U1,     LTI, Qc, Rinv,     tau, info)
            ! Compute Gamma = 0.5*(T0 @ (U1.T @ U0) + Tt @ (U1.T @ Ut))
            call mat_zero(QU); Swrk0 = 0.0_wp                           ! we use QU as a scratch array
            call mat_mult(Swrk0, X%U, U0); call mat_mult(QU, T0, Swrk0)
            call mat_mult(Swrk0, X%U, Ut); call mat_mult(T0, Tt, Swrk0) ! overwrite T0 with Gamma
            call mat_axpby(T0, 0.5_wp, QU, 0.5_wp)
            ! Update X to most recent value
            call mat_copy(X%U, U1)
            ! reverse steps based on Gamma
            call S_step_riccati(X, U1, QU, LTI, Qc, Rinv, 0.5*tau, info, reverse=.true., NL=T0)
            call K_step_riccati(X, U1, QU, LTI, Qc, Rinv, 0.5*tau, info, reverse=.true., NL=T0)
         end if
      else
         ! first order in time
         call K_step_riccati(   X, U1, QU, LTI, Qc, Rinv,     tau, info)
         call S_step_riccati(   X, U1, QU, LTI, Qc, Rinv,     tau, info)
         call L_step_riccati(   X, U1,     LTI, Qc, Rinv,     tau, info)
      end if
      
      !> Copy data to output
      call mat_copy(X%U, U1)
               
      return
   end subroutine G_forward_map_riccati

   subroutine K_step_riccati(X, U1, QU, LTI, Qc, Rinv, tau, info, reverse, NL)
      !> Low-rank state
      class(abstract_sym_low_rank_state), intent(inout) :: X         ! current state
      class(abstract_vector),             intent(out)   :: U1(:)     ! updated basis
      class(abstract_vector),             intent(inout) :: QU(:)     ! tmp
      ! LTI system
      class(abstract_lti_system),         intent(in)    :: LTI
      real(kind=wp),                      intent(in)    :: Qc(:,:)
      real(kind=wp),                      intent(in)    :: Rinv(:,:)
      !> Integration step size
      real(kind=wp),                      intent(in)    :: tau
      !> Information flag
      integer,                            intent(out)   :: info
      !> in second order scheme, are we in the forward or reverse branch?
      logical, optional,                  intent(in)    :: reverse
      logical                                           :: reverse_order
      !> precomputed non-linear term
      class(abstract_vector), optional,   intent(in)    :: NL(:)

      !> Local variables
      integer,                            allocatable   :: perm(:)   ! Permutation vector
      integer                                           :: rk

      ! Optional arguments
      reverse_order = optval(reverse, .false.)

      info = 0

      rk = size(X%U)
      if (.not. allocated(Uwrk0)) allocate(Uwrk0(1:rk), source=X%U(1)); 
      if (.not. allocated(Swrk0)) allocate(Swrk0(1:rk,1:rk)); 
      allocate(perm(1:rk)); perm = 0
      call mat_zero(Uwrk0); Swrk0 = 0.0_wp

      ! Constant part --> QU
      if (.not. reverse_order) then
         ! compute QU and pass to S step
         call apply_outerproduct_w(QU, X%U, LTI%CT, Qc)
      end if

      ! non-linear part --> Uwrk0
      call mat_mult(U1, X%U, X%S)                                  ! K0 = U0 @ S0
      if (.not.present(NL)) then
         call apply_p_outerproduct_w(Swrk0, X%U, U1, LTI%B, Rinv)  ! (U0.T) @ B @ R^(-1) @ B.T @ K0
         call mat_mult(Uwrk0, U1, Swrk0)                           ! K0 @ Swrk0
      else  ! non-linear term precomputed
         call mat_copy(Uwrk0, NL)
      end if
      
      ! combine to form G( K @ U.T ) @ U --> Uwrk0
      call mat_axpby(Uwrk0, -1.0_wp, QU, 1.0_wp)

      !> Construct solution U1
      call mat_axpby(U1, 1.0_wp, Uwrk0, tau)                       ! K0 + tau*Kdot

      !> Orthonormalize in-place
      call qr_factorization(U1, Swrk0, perm, info, ifpivot = .true.)
      call apply_permutation(Swrk0, perm, trans = .true.)
      X%S = Swrk0

      return
   end subroutine K_step_riccati

   subroutine S_step_riccati(X, U1, QU, LTI, Qc, Rinv, tau, info, reverse, NL)
      !> Low-rank state
      class(abstract_sym_low_rank_state), intent(inout) :: X      ! current state
      class(abstract_vector),             intent(in)    :: U1(:)  ! updated basis
      class(abstract_vector),             intent(inout) :: QU(:)  ! basis
      ! LTI system
      class(abstract_lti_system),         intent(in)    :: LTI
      real(kind=wp),                      intent(in)    :: Qc(:,:)
      real(kind=wp),                      intent(in)    :: Rinv(:,:)
      !> Integration step size
      real(kind=wp),                      intent(in)    :: tau
      !> Information flag
      integer,                            intent(out)   :: info
      !> in second order scheme, are we in the forward or reverse branch?
      logical, optional,                  intent(in)    :: reverse
      logical                                           :: reverse_order
      !> precomputed data
      class(abstract_vector), optional,   intent(in)    :: NL(:)

      !> Local variables
      integer                                           :: rk

      info = 0

      ! Optional arguments
      reverse_order = optval(reverse, .false.)

      rk = size(X%U)
      if (.not. allocated(Swrk0)) allocate(Swrk0(1:rk,1:rk))
      if (.not. allocated(Swrk1)) allocate(Swrk1(1:rk,1:rk))
      Swrk0 = 0.0_wp; Swrk1 = 0.0_wp

      ! Constant part --> Swrk0
      if (reverse_order) then
         ! compute QU and pass to K step
         call apply_outerproduct_w(QU, X%U, LTI%CT, Qc)
      endif
      call mat_mult(Swrk0, U1, QU)

      ! non-linear part --> Swrk1
      if (.not.present(NL)) then
         call apply_p_outerproduct_w(Swrk1, X%U, U1, LTI%B, Rinv)    !       U0.T @ B @ R^(-1) @ B.T @ U1
         Swrk1 = matmul(X%S, matmul(Swrk1, X%S))                     ! S0 @ (U0.T @ B @ R^(-1) @ B.T @ U1) @ S0
      else ! non-linear term precomputed
         call mat_mult(Swrk1, U1, NL)
      end if

      ! combine to form -U1.T @ G( U1 @ S @ U0.T ) @ U0
      call mat_axpby(Swrk0, -1.0_wp, Swrk1, 1.0_wp)

      !> Construct solution S1
      call mat_axpby(X%S, 1.0_wp, Swrk0, tau)

      return
   end subroutine S_step_riccati

   subroutine L_step_riccati(X, U1, LTI, Qc, Rinv, tau, info)
      !> Low-rank state
      class(abstract_sym_low_rank_state), intent(inout) :: X
      class(abstract_vector),             intent(in)    :: U1(:)   ! basis after Kstep
      ! LTI system
      class(abstract_lti_system),         intent(in)    :: LTI
      real(kind=wp),                      intent(in)    :: Qc(:,:)
      real(kind=wp),                      intent(in)    :: Rinv(:,:)
      !> Integration step size
      real(kind=wp),                      intent(in)    :: tau
      !> Information flag
      integer,                            intent(out)   :: info

      !> Local variables
      integer                                           :: rk

      info = 0

      rk = size(X%U)
      if (.not. allocated(Uwrk0)) allocate(Uwrk0(1:rk), source=X%U(1))
      if (.not. allocated(Uwrk1)) allocate(Uwrk1(1:rk), source=X%U(1))
      if (.not. allocated(Swrk0)) allocate(Swrk0(1:rk,1:rk)) 
      call mat_zero(Uwrk0); call mat_zero(Uwrk1); Swrk0 = 0.0_wp

      call mat_mult(Uwrk1, X%U, transpose(X%S))                   ! L0.T                                    U0 @ S.T

      ! Constant part --> Uwrk0
      call apply_outerproduct_w(Uwrk0, U1, LTI%CT, Qc)

      ! non-linear part --> U
      call apply_p_outerproduct_w(Swrk0, U1, Uwrk1, LTI%B, Rinv)  !               U1.T @ B @ R^(-1) @ B.T @ U0 @ S.T
      call mat_mult(X%U, Uwrk1, Swrk0)                            ! (U0 @ S.T) @ (U1.T @ B @ R^(-1) @ B.T @ U0 @ S.T)

      ! combine to form U1.T @ G( U1.T@L.T )
      call mat_axpby(Uwrk0, 1.0_wp, X%U, -1.0_wp)

      !> Construct solution U1
      call mat_axpby(Uwrk1, 1.0_wp, Uwrk0, tau)                   ! L0 + tau*Ldot

      !> Update S
      call mat_mult(X%S, Uwrk1, U1)

      return
   end subroutine L_step_riccati

end module lightROM_RiccatiSolvers
