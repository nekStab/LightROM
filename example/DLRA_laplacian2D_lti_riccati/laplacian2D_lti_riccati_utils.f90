module Laplacian2D_LTI_Riccati_Utils
   ! Standard Library.
   use stdlib_math, only : linspace
   use stdlib_optval, only : optval
   use stdlib_linalg, only : diag
   ! RKLIB module for time integration.
   use rklib_module
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov, only : wp => dp
   use LightKrylov_AbstractVectors ! zero_basis
   use LightKrylov_Utils, only : assert_shape
   ! Laplacian
   use Laplacian2D_LTI_Riccati_Base
   use laplacian2D_LTI_Riccati_Operators
   implicit none

   private
   ! initialisation
   public :: initialize_problem
   ! utils for state vector/matrix
   public :: get_state, set_state, init_rand
   ! initial condition
   public :: generate_random_initial_condition
   ! misc
   public :: CARE, sval
   
contains

   !---------------------------------------
   !-----     CONSTRUCT THE MESH      -----
   !---------------------------------------

   subroutine initialize_problem(magQ, magR)
      implicit none
      real(wp), intent(in)  :: magQ, magR
      ! internals
      real(wp), allocatable :: x(:)
      real(wp)              :: Bwrk(N,rk_b)
      integer :: i

      !> Construct mesh.
      x = linspace(-L/2, L/2, nx)

      ! Define C, Qc & compute CTQcC
      Qc = magQ*eye(rk_c)

      call init_rand(CT, ifnorm = .false.)
      call get_state(CTdata, CT)
      CTQcCdata = matmul(CTdata, matmul( Qc, transpose(CTdata)))
      CTQcC(1:N**2) = reshape(CTQcCdata, shape(CTQcC))

      ! Define B, Rinv & compule BRinvBT
      if (magR .lt. atol_dp) then
         Rinv = 0.0_wp
      else
         Rinv = 1/magR*eye(rk_b)
      endif

      call init_rand(B, ifnorm = .false.)
      Bdata = 0.0_wp
      call get_state(Bdata, B)
      Bwrk = 0.0_wp
      Bwrk = matmul(Bdata, Rinv)
      BRinvBTdata = matmul( Bwrk, transpose(Bdata) )

      return
   end subroutine initialize_problem

   !--------------------------------------------------------------------
   !-----     UTILITIES FOR STATE_VECTOR AND STATE MATRIX TYPES    -----
   !--------------------------------------------------------------------

   subroutine get_state(mat_out, state_in)
      !! Utility function to transfer data from a state vector to a real array
      real(wp),                   intent(out) :: mat_out(:,:)
      class(abstract_vector_rdp), intent(in)  :: state_in(:)
      ! internal variables
      integer :: k, kdim
      mat_out = 0.0_wp
      select type (state_in)
      type is (state_vector)
         kdim = size(state_in)
         call assert_shape(mat_out, (/ N, kdim /), 'get_state -> state_vector', 'mat_out')
         do k = 1, kdim
            mat_out(:,k) = state_in(k)%state
         end do
      type is (state_matrix)
         call assert_shape(mat_out, (/ N, N /), 'get_state -> state_matrix', 'mat_out')
         mat_out = reshape(state_in(1)%state, (/ N, N /))
      end select
      return
   end subroutine get_state

   subroutine set_state(state_out, mat_in)
      !! Utility function to transfer data from a real array to a state vector
      class(abstract_vector_rdp), intent(out) :: state_out(:)
      real(wp),                   intent(in)  :: mat_in(:,:)
      ! internal variables
      integer       :: k, kdim
      select type (state_out)
      type is (state_vector)
         kdim = size(state_out)
         call assert_shape(mat_in, (/ N, kdim /), 'set_state -> state_vector', 'mat_in')
         call zero_basis(state_out)
         do k = 1, kdim
            state_out(k)%state = mat_in(:,k)
         end do
      type is (state_matrix)
         call assert_shape(mat_in, (/ N, N /), 'set_state -> state_matrix', 'mat_in')
         call zero_basis(state_out)
         state_out(1)%state = reshape(mat_in, shape(state_out(1)%state))
      end select
      return
   end subroutine set_state

   subroutine init_rand(state, ifnorm)
      !! Utility function to initialize a state vector with random data
      class(abstract_vector_rdp), intent(inout)  :: state(:)
      logical, optional,          intent(in)     :: ifnorm
      ! internal variables
      integer :: k, kdim
      logical :: normalize
      normalize = optval(ifnorm,.true.)
      select type (state)
      type is (state_vector)
         kdim = size(state)
         do k = 1, kdim
            call state(k)%rand(ifnorm = normalize)
         end do
      type is (state_matrix)
         kdim = size(state)
         do k = 1, kdim
            call state(k)%rand(ifnorm = normalize)
         end do
      end select
      return
   end subroutine init_rand

   !--------------------------------------
   !-----     INITIAL CONDITIONS     -----
   !--------------------------------------

   subroutine generate_random_initial_condition(U, S, rk)
      class(state_vector),   intent(out) :: U(:)
      real(wp),              intent(out) :: S(:,:)
      integer,               intent(in)  :: rk
      ! internals
      class(state_vector),   allocatable :: Utmp(:)
      integer,               allocatable :: perm(:)
      ! SVD
      real(wp)                           :: U_svd(rk,rk)
      real(wp)                           :: S_svd(rk)
      real(wp)                           :: V_svd(rk,rk)
      integer                            :: i, info

      if (size(U) < rk) then
         write(*,*) 'Input krylov basis size incompatible with requested rank', rk
         STOP 1
      else
         call init_rand(U, .false.)
      end if
      if (size(S,1) < rk) then
         write(*,*) 'Input coefficient matrix size incompatible with requested rank', rk
         STOP 1
      else if (size(S,1) /= size(S,2)) then
         write(*,*) 'Input coefficient matrix must be square.'
         STOP 2
      else
         S = 0.0_wp
      end if
      ! perform QR
      allocate(perm(1:rk)); perm = 0
      allocate(Utmp(1:rk), source=U(1:rk))
      call qr(Utmp, S, perm, info, verbosity=.false.)
      if (info /= 0) write(*,*) '  [generate_random_initial_condition] Info: Colinear vectors detected in QR, column ', info
      ! perform SVD
      call svd(S(:,1:rk), U_svd(:,1:rk), S_svd(1:rk), V_svd(1:rk,1:rk))
      S = diag(S_svd)
      block
         class(abstract_vector_rdp), allocatable :: Xwrk(:)
         call linear_combination(Xwrk, Utmp, U_svd)
         call copy_basis(U, Xwrk)
      end block
      
   end subroutine

   !------------------------
   !-----     MISC     -----
   !------------------------

   function CARE(X,A,Q,BRinvBT) result(Y)
      real(wp), dimension(n,n) :: X, A, Q, BRinvBT, Y
      Y = matmul(transpose(A), X) + matmul(X, A) + Q - matmul(X, matmul(BRinvBT, X))
   end function CARE

   subroutine sval(X, svals)
      real(wp), intent(in) :: X(:,:)
      real(wp)             :: svals(min(size(X, 1), size(X, 2)))
      ! internals
      real(wp)             :: U(size(X, 1), size(X, 1))
      real(wp)             :: VT(size(X, 2), size(X, 2))
  
      ! Perform SVD
      call svd(X, U, svals, VT)
    
   end subroutine

end module Laplacian2D_LTI_Riccati_Utils