module LightROM_LyapunovSolvers
   use LightKrylov
   use LightKrylov_expmlib
   use LightKrylov_BaseKrylov
   use LightROM_AbstractLTIsystems
   use LightROM_LyapunovUtils
   use LightROM_utils
   use stdlib_linalg, only : eye
   use stdlib_optval, only : optval
   implicit none

   !> work arrays
   class(abstract_vector) , allocatable   :: U1(:)
   class(abstract_vector),  allocatable   :: Uwrk(:)
   class(abstract_vector) , allocatable   :: BBTU(:)
   real(kind=wp),           allocatable   :: Swrk(:,:)

   private
   public :: numerical_low_rank_splitting_lyapunov_integrator
   public :: M_forward_map, G_forward_map_lyapunov, K_step_lyapunov, S_step_lyapunov, L_step_lyapunov

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
   ! The algorithm constitutes a numerical integrator for the matrix-valued differential Lyapunov
   ! equation of the form
   !
   !          dX/dt = A @ X + X @ A.T + B.T @ B .
   !
   ! where A is a nxn Hurwitz matrix, X is SPD and B.T @ B is a rank-m rhs (m<<n).
   !
   ! Since A is Hurwitz, the equations converges to steady state for t -> infty, which corresponds
   ! to the associated algebrait Lyapunov equation of the form
   !
   !              0 = A @ X + X @ A.T + B.T @ B .
   !
   ! The algorithm is based on three ideas:
   ! a) The operator splitting scheme proposed by Lubich & Oseledets (2014) that splits the 
   !    right-hand side of the differential equation into a linear stiff part that is solved
   !    explicitly and a possibly non-linear non-stiff part which is solved numerically. The
   !    two operators are then composed to obtain the integrator for the full Lyapunov equation.
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
   ! - Separate integration of the stiff inhomogeneous part of the Lyapunov equation and the
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
   ! - The current implementation can only deal with constant linear inhomogeneities. The
   !   algorithm can be extended to deal with general (non-linear) inhomogeneities (see
   !   Lubich & Oseledets (2014) and Mena et al. (2018))
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
   ! - B      : Low-rank rhs of the Lyapunov equation     [Input]
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
   subroutine numerical_low_rank_splitting_lyapunov_integrator(X,LTI,Tend,tau,torder,info,exptA,iftrans)
      !> Low-rank state
      class(abstract_sym_low_rank_state),  intent(inout) :: X
      ! LTI system
      class(abstract_lti_system),          intent(in)    :: LTI
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
      integer                                :: istep, nsteps, iostep
      real(kind=wp)                          :: T
      logical, parameter                     :: verbose = .false.
      procedure(abstract_exptA), pointer     :: p_exptA => null()

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
         call numerical_low_rank_splitting_lyapunov_step(X, LTI, tau, torder, info, p_exptA, trans)

         T = T + tau
         !> here we should do some checks such as whether we have reached steady state
         if ( mod(istep,iostep) .eq. 0 ) then
            if (verbose) then
               write(*, *) "INFO : ", ISTEP, " steps of DLRA computed. T = ",T
            endif
         endif
      enddo dlra
      return

   end subroutine numerical_low_rank_splitting_lyapunov_integrator

   !-----------------------------
   !-----     UTILITIES     -----
   !-----------------------------
   subroutine numerical_low_rank_splitting_lyapunov_step(X, LTI, tau, torder, info, exptA, iftrans)
      !> Low-rank state
      class(abstract_sym_low_rank_state), intent(inout) :: X
      ! LTI system
      class(abstract_lti_system),         intent(in)    :: LTI
      !> Integration step size
      real(kind=wp),                      intent(in)    :: tau
      !> Order of time integration
      integer,                            intent(in)    :: torder
      !> Information flag
      integer,                            intent(out)   :: info
      !> Routine for computation of the exponential propagator
      procedure(abstract_exptA)                         :: exptA
      !> use transpose
      logical,                  optional, intent(in)    :: iftrans
      logical                                           :: trans

      !> Local variables
      integer        ::   istep, nsteps, integrator

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
         call M_forward_map(         X, LTI, tau, info, exptA, trans)
         call G_forward_map_lyapunov(X, LTI, tau, info)
      case (2) ! Strang splitting
         call M_forward_map(         X, LTI, 0.5*tau, info, exptA, trans)
         call G_forward_map_lyapunov(X, LTI,     tau, info)
         call M_forward_map(         X, LTI, 0.5*tau, info, exptA, trans)
      end select

      return

   end subroutine numerical_low_rank_splitting_lyapunov_step

   subroutine M_forward_map(X, LTI, tau, info, exptA, iftrans)
      !> Low-rank state
      class(abstract_sym_low_rank_state), intent(inout) :: X
      ! LTI system
      class(abstract_lti_system),         intent(in)    :: LTI
      !> Integration step size
      real(kind=wp),                      intent(in)    :: tau
      !> Information flag
      integer,                            intent(out)   :: info
      !> Routine for computation of the exponential propagator
      procedure(abstract_exptA)                         :: exptA
      !> use transpose
      logical, optional,                  intent(in)    :: iftrans
      logical                                           :: trans

      !> Local variables
      class(abstract_vector),             allocatable   :: Uwrk  ! basis
      real(kind=wp),                      allocatable   :: R(:,:)   ! QR coefficient matrix
      integer,                            allocatable   :: perm(:)   ! Permutation vector
      real(kind=wp),                      allocatable   :: wrk(:,:)
      integer                                           :: i, rk

      ! Optional argument
      trans = optval(iftrans, .false.)

      rk = size(X%U)
      !> Allocate memory
      allocate(R(1:rk,1:rk));   R = 0.0_wp 
      allocate(perm(1:rk));     perm = 0
      allocate(wrk(1:rk,1:rk)); wrk = 0.0_wp

      !> Apply propagator to initial basis
      if (.not. allocated(Uwrk)) allocate(Uwrk, source=X%U(1))
      call Uwrk%zero()
      do i = 1, rk
         call exptA(Uwrk, LTI%A, X%U(i), tau, info, trans)
         call X%U(i)%axpby(0.0_wp, Uwrk, 1.0_wp) ! overwrite old solution
      enddo
      !> Reorthonormalize in-place
      call qr_factorization(X%U, R, perm, info, ifpivot = .true.)
      !> Update low-rank coefficient matrix
      call apply_permutation(R, perm, trans = .true.)
      wrk = matmul(X%S, transpose(R))
      X%S = matmul(R, wrk)

      return
   end subroutine M_forward_map

   subroutine G_forward_map_lyapunov(X, LTI, tau, info)
      !> Low-rank state
      class(abstract_sym_low_rank_state), intent(inout) :: X
      ! LTI system
      class(abstract_lti_system),         intent(in)    :: LTI
      !> Integration step size
      real(kind=wp),                      intent(in)    :: tau
      !> Information flag
      integer,                            intent(out)   :: info

      !> Local variables
      class(abstract_vector),             allocatable   :: U1(:)
      class(abstract_vector),             allocatable   :: BBTU(:)
      integer                                           :: rk

      rk = size(X%U)
      if (.not. allocated(U1))   allocate(U1(1:rk),   source=X%U(1))
      if (.not. allocated(BBTU)) allocate(BBTU(1:rk), source=X%U(1))
      call mat_zero(U1); call mat_zero(BBTU)

      call K_step_lyapunov(X, U1, BBTU, LTI%B, tau, info)

      call S_step_lyapunov(X, U1, BBTU,        tau, info)

      call L_step_lyapunov(X, U1,       LTI%B, tau, info)
      
      !> Copy data to output
      call mat_copy(X%U, U1)
               
      return
   end subroutine G_forward_map_lyapunov

   subroutine K_step_lyapunov(X, U1, BBTU, B, tau, info)
      !> Low-rank state
      class(abstract_sym_low_rank_state), intent(inout) :: X
      class(abstract_vector),             intent(out)   :: U1(:)  ! basis
      class(abstract_vector),             intent(out)   :: BBTU(:)  ! basis
      !> low-rank rhs
      class(abstract_vector),             intent(in)    :: B(:)
      !> Integration step size
      real(kind=wp),                      intent(in)    :: tau
      !> Information flag
      integer,                            intent(out)   :: info

      !> Local variables
      real(kind=wp),                      allocatable   :: Swrk(:,:)
      integer,                            allocatable   :: perm(:)   ! Permutation vector
      integer                                           :: rk

      info = 0

      rk = size(X%U)
      if (.not. allocated(Swrk)) allocate(Swrk(1:rk,1:rk));
      Swrk = 0.0_wp
      allocate(perm(1:rk)); perm = 0

      call mat_mult(U1, X%U, X%S)             ! K0
      call apply_outerproduct(BBTU, B, X%U)   ! Kdot
      !> Construct solution U1
      call mat_axpby(U1, 1.0_wp, BBTU, tau) ! K0 + tau*Kdot
      !> Orthonormalize in-place
      call qr_factorization(U1, Swrk, perm, info, ifpivot = .true.)
      call apply_permutation(Swrk, perm, trans = .true.)
      X%S = Swrk

      return
   end subroutine K_step_lyapunov

   subroutine S_step_lyapunov(X, U1, BBTU, tau, info)
      !> Low-rank state
      class(abstract_sym_low_rank_state), intent(inout) :: X
      class(abstract_vector),             intent(in)    :: U1(:)    ! updated basis
      class(abstract_vector),             intent(in)    :: BBTU(:)  ! basis
      !> Integration step size
      real(kind=wp),                      intent(in)    :: tau
      !> Information flag
      integer,                            intent(out)   :: info

      !> Local variables
      real(kind=wp),                      allocatable   :: Swrk(:,:)
      integer                                           :: rk

      info = 0

      rk = size(X%U)
      if (.not. allocated(Swrk)) allocate(Swrk(1:rk,1:rk))
      Swrk = 0.0_wp
      call mat_mult(Swrk, U1, BBTU)          ! - Sdot
      !> Construct solution S1
      call mat_axpby(X%S, 1.0_wp, Swrk, -tau)

      return
   end subroutine S_step_lyapunov

   subroutine L_step_lyapunov(X, U1, B, tau, info)
      !> Low-rank state
      class(abstract_sym_low_rank_state), intent(inout) :: X
      class(abstract_vector),             intent(in)    :: U1(:)   ! basis after Kstep
      !> low-rank rhs
      class(abstract_vector),             intent(in)    :: B(:)
      !> Integration step size
      real(kind=wp),                      intent(in)    :: tau
      !> Information flag
      integer,                            intent(out)   :: info

      !> Local variables
      class(abstract_vector),             allocatable   :: Uwrk(:)
      integer                                           :: rk

      info = 0

      rk = size(X%U)
      if (.not. allocated(Uwrk)) allocate(Uwrk(1:rk), source=X%U(1))
      call mat_zero(Uwrk)

      call mat_mult(Uwrk, X%U, transpose(X%S))  ! L0.T
      call apply_outerproduct(X%U, B, U1)     ! Ldot.T
      !> Construct solution Uwrk.T
      call mat_axpby(Uwrk, 1.0_wp, X%U, tau)
      !> Update S
      call mat_mult(X%S, Uwrk, U1)

      return
   end subroutine L_step_lyapunov

end module lightROM_LyapunovSolvers
