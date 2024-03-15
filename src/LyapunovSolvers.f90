module LightROM_LyapunovSolvers
   use LightKrylov
   use LightROM_LyapunovUtils
   use LightROM_utils
   use LightROM_expmlib
   use stdlib_linalg, only : eye
   implicit none

   !> work arrays
   class(abstract_vector),  allocatable   :: Uwrk(:)
   real(kind=wp),           allocatable   :: Swrk(:,:)

   private
   public :: numerical_low_rank_splitting_integrator, M_forward_map, G_forward_map

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
   ! - The current implementation is first order in time for the non-stiff solution. 
   !   Higher order schemes are possible (see Lubich & Oseledets (2014))
   ! - The current implementation does not require an adjoint integrator. This means that
   !   the temporal order of the basic operator splitting scheme is limited to 1 (Lie-Trotter
   !   splitting) or at most 2 (Strang splitting). Higher order integrators are possible, but 
   !   require at least some backward integration (via the adjoint) in BOTH parts. 
   !   (see Sheng-Suzuki and Goldman-Kaper theorems)
   !
   ! Input/Output Parameters:
   ! ------------------------
   ! - A        : Linear Operator (SPD)                     [Input]
   ! - B        : Low-rank rhs of the Lyapunov equation     [Input]
   ! - X        : Initial/Updated solution                  [Input/Output]
   ! - info     : Iteration Information flag                [Output]
   ! - tol      : Tolerance for convergence                 [Optional, Input]
   ! - verbosity: Verbosity control flag                    [Optional, Input]
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
   subroutine numerical_low_rank_splitting_integrator(U,S,A,B,Tend,tau,torder,info)
      !> Low-rank factors
      class(abstract_vector) , intent(inout) :: U(:) ! basis
      real(kind=wp)          , intent(inout) :: S(:,:) ! coefficients
      !> System
      class(abstract_linop)  , intent(in)    :: A    ! linear operator
      class(abstract_vector) , intent(in)    :: B(:) !
      !> Integration time and step size
      real(kind=wp)          , intent(in)    :: Tend
      real(kind=wp)          , intent(inout) :: tau  ! desired/effective time step
      !> Order of time integration
      integer                , intent(in)    :: torder
      !> Information flag
      integer                , intent(out)   :: info

      !> Local variables   
      integer                                :: rk, istep, nsteps, iostep
      real(kind=wp)                          :: T
      logical, parameter                     :: verbose = .false.
      T = 0.0_wp
      rk = size(U)

      ! --> Compute number of steps
      nsteps = floor(Tend/tau)

      iostep = nsteps/10
      if ( iostep .eq. 0 ) then
         iostep = 10
      endif

      dlra : do istep = 1, nsteps
         !write(*,*) 'istep', istep
         !> dynamical low-rank approximation step
         call numerical_low_rank_splitting_step(U, S, A, B, tau, torder, info)

         T = T + tau
         !> here we should do some checks such as whether we have reached steady state
         if ( mod(istep,iostep) .eq. 0 ) then
            if (verbose) then
               write(*, *) "INFO : ", ISTEP, " steps of DLRA computed. T = ",T
            endif
         endif
      enddo dlra
      return

      end subroutine numerical_low_rank_splitting_integrator

   !contains

      !-----------------------------
      !-----     UTILITIES     -----
      !-----------------------------
      subroutine numerical_low_rank_splitting_step(U, S, A, B, tau, torder, info)
         !> Low-rank factors
         class(abstract_vector) , intent(inout) :: U(:) ! basis
         real(kind=wp)          , intent(inout) :: S(:,:) ! coefficients
         !> Linear operator and rhs Cholesky factor
         class(abstract_linop)  , intent(in)    :: A    ! linear operator
         class(abstract_vector) , intent(in)    :: B(:)
         !> Integration step size
         real(kind=wp)          , intent(in)    :: tau
         !> Order of time integration
         integer                , intent(in)    :: torder
         !> Information flag
         integer                , intent(out)   :: info

         !> Local variables
         integer        ::   istep, nsteps, integrator

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
            call M_forward_map(U, S, A, tau, info)
            call G_forward_map(U, S, B, tau, info)
         case (2) ! Strang splitting
            call M_forward_map(U, S, A, 0.5*tau, info)
            call G_forward_map(U, S, B,     tau, info)
            call M_forward_map(U, S, A, 0.5*tau, info)
         end select

         return

      end subroutine numerical_low_rank_splitting_step

      subroutine M_forward_map(U, S, A, tau, info)
         !> Low-rank factors
         class(abstract_vector) , intent(inout) :: U(:)   ! basis
         real(kind=wp)          , intent(inout) :: S(:,:) ! coefficients
         !> Linear operator
         class(abstract_linop)  , intent(in)    :: A 
         !> Integration step size
         real(kind=wp)          , intent(in)    :: tau
         !> Information flag
         integer                , intent(out)   :: info

         !> Local variables
         class(abstract_vector) , allocatable   :: Uwrk  ! basis
         real(kind=wp)          , allocatable   :: R(:,:)   ! QR coefficient matrix
         real(kind=wp)          , allocatable   :: wrk(:,:)
         integer                                :: i, rk

         rk = size(U)
         !> Allocate memory
         allocate(R(1:rk,1:rk)); allocate(wrk(1:rk,1:rk)); 
         R = 0.0_wp; wrk = 0.0_wp
         call print_mat(rk,rk,S,'M step start')

         !> Apply propagator to initial basis
         if (.not. allocated(Uwrk)) allocate(Uwrk, source=U(1))
         call Uwrk%zero()
         do i = 1, rk
            call A%matvec(U(i), Uwrk)
            call U(i)%axpby(0.0_wp, Uwrk, 1.0_wp) ! overwrite old solution
         enddo
         call print_mat(rk,rk,S,'M step pr orth')
         !> Reorthonormalize in-place
         call qr_factorization(U, R, info)
         !> Update low-rank coefficient matrix
         wrk = matmul(S, transpose(R))
         S   = matmul(R, wrk)

         call print_mat(rk,rk,S,'M step')

         return
      end subroutine M_forward_map

      subroutine G_forward_map(U, S, B, tau, info)
         !> Low-rank factors
         class(abstract_vector) , intent(inout) :: U(:)   ! basis
         real(kind=wp)          , intent(inout) :: S(:,:) ! coefficients
         !> low-rank rhs
         class(abstract_vector) , intent(in)    :: B(:)
         !> Integration step size
         real(kind=wp)          , intent(in)    :: tau
         !> Information flag
         integer                , intent(out)   :: info

         !> Local variables
         class(abstract_vector) , allocatable   :: U1(:) 
         integer                                :: rk

         rk = size(U)
         allocate(U1(1:rk), source=U(1)); call mat_zero(U1)

         call K_step(U1, S, U,     B, tau, info)

         call S_step(    S, U, U1, B, tau, info)

         call L_step(    S, U, U1, B, tau, info)
         
         !> Copy data to output
         call mat_copy(U, U1)
                  
         return
      end subroutine G_forward_map

      subroutine K_step(U1, S, U, B, tau, info)
         !> Low-rank factors
         class(abstract_vector) , intent(out)   :: U1(:)  ! basis
         real(kind=wp)          , intent(inout) :: S(:,:) ! coefficients
         !> Low-rank factors
         class(abstract_vector) , intent(in)    :: U(:)   ! basis
         
         !> low-rank rhs
         class(abstract_vector) , intent(in)    :: B(:)
         !> Integration step size
         real(kind=wp)          , intent(in)    :: tau
         !> Information flag
         integer                , intent(out)   :: info

         !> Local variables
         class(abstract_vector) , allocatable   :: Uwrk(:)   
         real(kind=wp)          , allocatable   :: Swrk(:,:)  
         integer                                :: rk

         info = 0

         rk = size(U)
         if (.not. allocated(Uwrk)) allocate(Uwrk(1:rk), source=U(1));
         if (.not. allocated(Swrk)) allocate(Swrk(1:rk,1:rk)); 
         call mat_zero(Uwrk); Swrk = 0.0_wp

         call mat_mult(U1, U, S)               ! K0
         call apply_outerproduct(Uwrk, B, U)   ! Kdot
         !> Construct solution U1
         call mat_axpby(U1, 1.0_wp, Uwrk, tau) ! K0 + tau*Kdot
         !> Orthonormalize in-place
         call qr_factorization(U1, Swrk, info)
         S = Swrk

         return
      end subroutine K_step

      subroutine S_step(S, U, U1, B, tau, info)
         !> Low-rank factors
         real(kind=wp)          , intent(inout) :: S(:,:) ! coefficients
         !> Low-rank factors
         class(abstract_vector) , intent(in)    :: U(:)   ! old basis
         class(abstract_vector) , intent(in)    :: U1(:)  ! updated basis
         !> low-rank rhs
         class(abstract_vector) , intent(in)    :: B(:)
         !> Integration step size
         real(kind=wp)          , intent(in)    :: tau
         !> Information flag
         integer                , intent(out)   :: info

         !> Local variables
         class(abstract_vector) , allocatable   :: Uwrk(:)
         real(kind=wp)          , allocatable   :: Swrk(:,:)
         integer                                :: rk

         info = 0

         rk = size(U)
         if (.not. allocated(Uwrk)) allocate(Uwrk(1:rk), source=U(1));
         if (.not. allocated(Swrk)) allocate(Swrk(1:rk,1:rk)); 
         call mat_zero(Uwrk); Swrk = 0.0_wp

         call apply_outerproduct(Uwrk, B, U)
         call mat_mult(Swrk, U1, Uwrk)          ! - Sdot
         !> Construct solution S1
         call mat_axpby(S, 1.0_wp, Swrk, -tau)

         return
      end subroutine S_step

      subroutine L_step(S, U, U1, B, tau, info)
         !> Low-rank factors
         real(kind=wp)          , intent(inout) :: S(:,:) ! coefficients
         !> Low-rank factors
         class(abstract_vector) , intent(inout) :: U(:)   ! basis before Kstep
         class(abstract_vector) , intent(in)    :: U1(:)   ! basis after Kstep
         !> low-rank rhs
         class(abstract_vector) , intent(in)    :: B(:)
         !> Integration step size
         real(kind=wp)          , intent(in)    :: tau
         !> Information flag
         integer                , intent(out)   :: info

         !> Local variables
         class(abstract_vector) , allocatable   :: Uwrk(:)
         integer                                :: rk

         info = 0

         rk = size(U)
         if (.not. allocated(Uwrk)) allocate(Uwrk(1:rk), source=U(1));
         call mat_zero(Uwrk)

         call mat_mult(Uwrk, U, transpose(S))  ! L0.T
         call apply_outerproduct(U, B, U1)     ! Ldot.T
         !> Construct solution Uwrk.T
         call mat_axpby(Uwrk, 1.0_wp, U, tau)
         !> Update S
         call mat_mult(S, Uwrk, U1)

         return
      end subroutine L_step 

   !end subroutine numerical_low_rank_splitting_integrator

end module lightROM_LyapunovSolvers
