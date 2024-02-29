module LightROM_LyapunovSolvers
   use LightKrylov
   use LightROM_LyapunovUtils
   !include "dtypes.h"

   private
   public :: numerical_low_rank_splitting_integrator

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
      integer :: istep, nsteps
      real(kind=wp) :: T
      T = 0.0_wp    

      ! --> Reset desired tau to match Tend
      nsteps = ceiling(Tend/tau)
      tau    = Tend/nsteps

      dlra : do istep = 1, nsteps
         write(*,*) 'DLRA',istep
         !> dynamical low-rank approximation step
         call numerical_low_rank_splitting_step(U, S, A, B, tau, torder, info)

         T = T + tau
         !> here we should do some checks such as whether we have reached steady state
         !if ( modulo(istep,iostep) .eq. 0 ) then
         !   if (verbose) then
         !      write(*, *) "INFO : ", IOSTEP, " steps of DLRA computed."
         !   endif
         !endif
      enddo dlra
      return

   contains

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
         write(*,*) '    step: integrator', integrator

         select case (integrator)
         case (1) ! Lie-Trotter splitting
            write(*,*) '    step: splitting: Lie-Trotter'
            call M_forward_map(U, S, A, tau, info)
            call G_forward_map(U, S, B, tau, info)
         case (2) ! Strang splitting
            call M_forward_map(U, S, A, 0.5*tau, info)
            call G_forward_map(U, S, B, tau, info)
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
         !> tolerance for the exponential propagator
         real(kind=wp), parameter :: tol = 1e-9

         !> Local variables
         class(abstract_vector) , allocatable   :: Uwrk(:)  ! basis
         real(kind=wp)          , allocatable   :: R(:,:)   ! QR coefficient matrix
         real(kind=wp)          , allocatable   :: wrk(:,:)
         integer                                :: i, rk

         rk = size(U)
         allocate(R(1:rk,1:rk)); R = 0.0_wp
         call mat_mult(R,U(1:rk),U(1:rk))
         write(*,*) '         U orthonormality (before M step):'
         do i = 1, rk
            write(*,'(4F8.3)') R(i,1:rk)
         enddo
         !> Apply propagator to initial basis
         allocate(Uwrk(1:rk), source=U(1)) ; call mat_zero(Uwrk)
         do i = 1, rk
            call A%matvec(U(i), Uwrk(i))
            !call U(i)%axpby(0.0_wp, Uwrk, 1.0_wp) ! overwrite old solution
         enddo
         call mat_mult(R,Uwrk(1:rk),Uwrk(1:rk))
         write(*,*) '         U orthonormality (after M step):'
         do i = 1, rk
            write(*,'(4E15.8)') R(i,1:rk)
         enddo
         !> Reorthonormalize in-place
         call qr_factorization(U, R, info)
         !> Update low-rank coefficient matrix
         write(*,*) '         R (after M step):'
         do i = 1, rk
            write(*,'(4F8.3)') R(i,1:rk)
         enddo
         allocate(wrk(1:rk,1:rk)); wrk = 0.0_wp
         wrk = matmul(S, transpose(R))
         S   = matmul(R, wrk)
         deallocate(Uwrk)
         deallocate(wrk)
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
         class(abstract_vector) , allocatable   :: Uwrk(:)
         real(kind=wp)          , allocatable   :: S1(:,:)
         real(kind=wp)          , allocatable   :: Swrk(:,:)
         integer :: i, rk

         rk = size(U)
         !> Allocate temporary arrays
         allocate(U1(1:rk), source=U(1)); allocate(Uwrk(1:rk), source=U(1))
         call mat_zero(U1); call mat_zero(Uwrk)
         allocate(Swrk(1:rk,1:rk)); allocate(S1(1:rk,1:rk)); 
         Swrk = 0.0_wp; S1 = 0.0_wp
         !> Step K equation:      
         !>       K0   = U @ S
         !>       Kdot = B @ B.T @ U
         call mat_mult(U1, U, S)     ! K0
         call apply_outerproduct(Uwrk, B, U)  ! Kdot
         !> Construct solution U1
         call mat_axpby(U, 1.0_wp, Uwrk, tau)
         !> Orthonormalize in-place
         call qr_factorization(U1, S1, info)
         !
         !> Step S equation
         !>       S0   = S1
         !>       Sdot = -U1.T @ B @ B.T @ UA
         call apply_outerproduct(Uwrk, B, U)
         call mat_mult(Swrk, U1, Uwrk)          ! - Sdot
         !> Construct solution S1
         call mat_axpby(S1, 1.0_wp, Swrk, -tau)
         !
         !> Step L equation:
         !>       L0   = St @ UA.T           or      L0.T   = UA @ St.T  
         !>       Ldot = U1.T @ B @ B.T      or      Ldot.T = B @ B.T @ U1
         call mat_mult(Uwrk, U, transpose(S1))  ! L0.T
         !> note: we no longer need the old basis U so we overwrite it
         call apply_outerproduct(U, B, U1)               ! Ldot.T
         !> Construct solution Uwrk.T
         call mat_axpby(Uwrk, 1.0_wp, U, tau)
         !> Update S
         call mat_mult(S, Uwrk, U1)
         !> Copy U1 --> U for output
         call mat_copy(U, U1)
         !> Deallocate work arrays
         deallocate(U1)
         deallocate(Uwrk)
         deallocate(S1)
         deallocate(Swrk)
         return
      end subroutine G_forward_map

      !subroutine K_step(K1,U,S,tau)
      !end subroutine K_step
!
      !subroutine S_step 
      !end subroutine S_step
!
      !subroutine L_step 
      !end subroutine L_step 

   end subroutine numerical_low_rank_splitting_integrator
end module lightROM_LyapunovSolvers
