module Laplacian2D_LTI_Lyapunov_Utils
   ! Standard Library.
   use stdlib_math, only : linspace
   use stdlib_optval, only : optval
   use stdlib_linalg, only : eye, diag, svd, svdvals
   ! RKLIB module for time integration.
   use rklib_module
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov_Constants, only : one_rdp, zero_rdp
   use LightKrylov_AbstractVectors ! linear_combination
   ! Laplacian
   use Laplacian2D_LTI_Lyapunov_Base
   use laplacian2D_LTI_Lyapunov_Operators
   implicit none

   private:: this_module
   ! mesh
   public :: initialize_mesh
   ! utilities for state matrix
   public :: get_state, set_state, init_rand
   ! initial conditions
   public :: generate_random_initial_condition
   ! misc
   public :: CALE, build_operator, reconstruct_TQ, solve_lyapunov

   character(len=*), parameter :: this_module = 'Laplacian2D_LTI_Lyapunov_Utils'

contains

   !---------------------------------------
   !-----     CONSTRUCT THE MESH      -----
   !---------------------------------------

   subroutine initialize_mesh()
      implicit none
      !> Mesh array.
      real(dp), allocatable :: x(:)
      integer :: i

      !> Construct mesh.
      x = linspace(-L/2, L/2, nx)

   end subroutine initialize_mesh

   !--------------------------------------------------------------------
   !-----     UTILITIES FOR STATE_VECTOR AND STATE MATRIX TYPES    -----
   !--------------------------------------------------------------------

   subroutine get_state(mat_out, state_in, procedure)
      !! Utility function to transfer data from a state vector to a real array
      real(dp),                   intent(out) :: mat_out(:,:)
      class(abstract_vector_rdp), intent(in)  :: state_in(:)
      character(len=*),           intent(in)  :: procedure
      ! internal variables
      character(len=*), parameter :: this_procedure = 'get_state'
      integer :: k, kdim
      mat_out = 0.0_dp
      select type (state_in)
      type is (state_vector)
         kdim = size(state_in)
         call assert_shape(mat_out, [ N, kdim ], 'mat_out', this_module, this_procedure//': '//trim(procedure))
         do k = 1, kdim
            mat_out(:,k) = state_in(k)%state
         end do
      type is (state_matrix)
         call assert_shape(mat_out, [ N, N ], 'mat_out', this_module, this_procedure//': '//trim(procedure))
         mat_out = reshape(state_in(1)%state, [ N, N ])
      class default
         call type_error('state', 'state_vector or state_matrix', 'IN', this_module, this_procedure)
      end select
   end subroutine get_state

   subroutine set_state(state_out, mat_in, procedure)
      !! Utility function to transfer data from a real array to a state vector
      class(abstract_vector_rdp), intent(out) :: state_out(:)
      real(dp),                   intent(in)  :: mat_in(:,:)
      character(len=*),           intent(in)  :: procedure
      ! internal variables
      character(len=*), parameter :: this_procedure = 'set_state'
      integer       :: k, kdim
      select type (state_out)
      type is (state_vector)
         kdim = size(state_out)
         call assert_shape(mat_in, [ N, kdim ], 'mat_in', this_module,this_procedure//': '//trim(procedure))
         call zero_basis(state_out)
         do k = 1, kdim
            state_out(k)%state = mat_in(:,k)
         end do
      type is (state_matrix)
         call assert_shape(mat_in, [ N, N ], 'mat_in', this_module, this_procedure//': '//trim(procedure))
         call zero_basis(state_out)
         state_out(1)%state = reshape(mat_in, shape(state_out(1)%state))
      class default
         call type_error('state', 'state_vector or state_matrix', 'IN', this_module, this_procedure)
      end select
   end subroutine set_state

   subroutine init_rand(state, ifnorm)
      !! Utility function to initialize a state vector with random data
      class(abstract_vector_rdp), intent(inout)  :: state(:)
      logical, optional,          intent(in)     :: ifnorm
      ! internal variables
      character(len=*), parameter :: this_procedure = 'init_rand'
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
      class default
         call type_error('state', 'state_vector or state_matrix', 'INOUT', this_module, this_procedure)
      end select
   end subroutine init_rand

   !--------------------------------------
   !-----     INITIAL CONDITIONS     -----
   !--------------------------------------

   subroutine generate_random_initial_condition(U, S, rk)
      class(state_vector),   intent(out) :: U(:)
      real(dp),              intent(out) :: S(:,:)
      integer,               intent(in)  :: rk
      ! internals
      character(len=*), parameter :: this_procedure = 'generate_random_initial_condition'
      class(state_vector),   allocatable :: Utmp(:)
      integer                            :: i, info
      real(dp) :: mean, std
      character(len=128) :: msg

      ! sanity check
      if (size(U) < rk) then
         write(msg,'(A,I0)') 'Input krylov basis size incompatible with requested rank ', rk
         call stop_error(msg, this_module, this_procedure)
      end if

      call zero_basis(U)
      call initialize_random_orthonormal_basis(U(:rk))
      call assert_shape(S, [ rk,rk ], 'S', this_module, this_procedure)
      S = zero_rdp
      mean = zero_rdp
      std  = one_rdp
      do i = 1, rk
         S(i,i) = normal(mean,std)
      end do
      write(msg,'(A,I0,A,I0,A)') 'size(U) = [ ', size(U),' ]: filling the first ', rk, ' columns with orthonormal noise.'
      call logger%log_information(msg, this_module, this_procedure)
   end subroutine

   !------------------------
   !-----     MISC     -----
   !------------------------

   function CALE(X,A,Q) result(res)
      real(dp), dimension(n,n) :: X, A, Q, res
      res = matmul(A, X) + matmul(X, transpose(A)) + Q
   end function CALE

   subroutine build_operator(A)
      !! Build the two-dimensional Laplace operator explicitly
      real(dp), intent(out) :: A(N,N)
      integer :: i, j, k

      A = -4.0_dp/dx2*eye(N)
      do i = 1, nx
         do j = 1, nx - 1
            k = (i-1)*nx + j
            A(k + 1, k) = 1.0_dp/dx2
            A(k, k + 1) = 1.0_dp/dx2
         end do 
      end do
      do i = 1, N-nx
         A(i, i + nx) = 1.0_dp/dx2
         A(i + nx, i) = 1.0_dp/dx2
      end do
   end subroutine build_operator

   subroutine reconstruct_TQ(T, Q, A, D, E, tw)
      !! Reconstruct tridiagonal matrix T and orthogonal projector Q from dsytd2 output (A, D, E)
      real(dp), intent(out) :: T(N,N)
      real(dp), intent(out) :: Q(N,N)
      real(dp), intent(in)  :: A(N,N)
      real(dp), intent(in)  :: D(N)
      real(dp), intent(in)  :: E(N-1)
      real(dp), intent(in)  :: tw(N-1)

      ! internal variables
      real(dp)  :: Hi(N,N)
      real(dp)  :: vec(N,1)
      integer :: i

      ! Build orthogonal Q = H(1) @  H(2) @ ... @ H(n-1)
      Q = eye(N)
      do i = 1, N - 1
         vec          = 0.0_dp
         vec(i+1,1)   = 1.0_dp
         vec(i+2:N,1) = A(i+2:N,i)
         Hi           = eye(N) - tw(i) * matmul( vec, transpose(vec) )
         Q            = matmul( Q, Hi )
      end do

      ! Build tridiagonal T
      T = 0.0_dp
      do i = 1, N
         T(i,i) = D(i)
      end do
      do i = 1, N - 1
         T(i,i+1) = E(i)
         T(i+1,i) = E(i)
      end do

   end subroutine reconstruct_TQ

   subroutine solve_lyapunov(X, A, P)
      !! Solve the Lyapunov equation directly
      real(dp), intent(out) :: X(N,N)
      !! Solution
      real(dp), intent(in)  :: A(N,N)
      !! Operator
      real(dp), intent(in)  :: P(N,N)
      !! Inhomogeneity
      ! Internal
      real(dp), dimension(N,N) :: T, Q, Z, V, W, Y
      real(dp), dimension(N)   :: Dm, wrk, wr, wi
      real(dp), dimension(N-1) :: E, tw
      real(dp)                 :: scale
      integer                  :: isgn, info

      ! Transform operator to tridiagonal form
      call dsytd2('L', N, A, N, Dm, E, tw, info)

      ! Reconstruct T and Q
      call reconstruct_TQ(T, Q, A, Dm, E, tw)

      ! compute real Schur form A = Z @ T @ Z.T
      call dhseqr('S', 'I', N, 1, N, T, N, wr, wi, Z, N, wrk, N, info )

      ! Change RHS Basis: base --> Q --> Z
      V = matmul(transpose(Q), matmul(-P, Q))
      W = matmul(transpose(Z), matmul( V, Z))

      ! Compute solution of Lyapunov equation for Schur decomposition
      isgn = 1; scale = 0.1_dp
      call dtrsyl('N', 'T', isgn, N, N, T, N, T, N, W, N, scale, info)

      ! Return to original basis to obtain X_ref: Z --> Q --> base
      Y = matmul(Z, matmul(W, transpose(Z)))
      X = matmul(Q, matmul(Y, transpose(Q)))

   end subroutine solve_lyapunov
   
end module Laplacian2D_LTI_Lyapunov_Utils
