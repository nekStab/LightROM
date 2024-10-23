module Laplacian2D_LTI_Lyapunov_Utils
   ! Standard Library.
   use stdlib_math, only : linspace
   use stdlib_optval, only : optval
   use stdlib_linalg, only : eye, diag, svd
   ! RKLIB module for time integration.
   use rklib_module
   ! LightKrylov for linear algebra.
   use LightKrylov, only : wp => dp
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
   public :: CALE, build_operator, reconstruct_TQ

   character(len=128), parameter :: this_module = 'Laplacian2D_LTI_Lyapunov_Utils'

contains

   !---------------------------------------
   !-----     CONSTRUCT THE MESH      -----
   !---------------------------------------

   subroutine initialize_mesh()
      implicit none
      !> Mesh array.
      real(wp), allocatable :: x(:)
      integer :: i

      !> Construct mesh.
      x = linspace(-L/2, L/2, nx)

      return
   end subroutine initialize_mesh

   !--------------------------------------------------------------------
   !-----     UTILITIES FOR STATE_VECTOR AND STATE MATRIX TYPES    -----
   !--------------------------------------------------------------------

   subroutine get_state(mat_out, state_in, procedure)
      !! Utility function to transfer data from a state vector to a real array
      real(wp),                   intent(out) :: mat_out(:,:)
      class(abstract_vector_rdp), intent(in)  :: state_in(:)
      character(len=*),           intent(in)  :: procedure
      ! internal variables
      integer :: k, kdim
      mat_out = 0.0_wp
      select type (state_in)
      type is (state_vector)
         kdim = size(state_in)
         call assert_shape(mat_out, [ N, kdim ], 'mat_out', this_module, 'get_state -> state_vector')
         do k = 1, kdim
            mat_out(:,k) = state_in(k)%state
         end do
      type is (state_matrix)
         call assert_shape(mat_out, [ N, N ], 'mat_out', this_module, 'get_state -> state_matrix')
         mat_out = reshape(state_in(1)%state, [ N, N ])
      end select
      return
   end subroutine get_state

   subroutine set_state(state_out, mat_in, procedure)
      !! Utility function to transfer data from a real array to a state vector
      class(abstract_vector_rdp), intent(out) :: state_out(:)
      real(wp),                   intent(in)  :: mat_in(:,:)
      character(len=*),           intent(in)  :: procedure
      ! internal variables
      integer       :: k, kdim
      select type (state_out)
      type is (state_vector)
         kdim = size(state_out)
         call assert_shape(mat_in, [ N, kdim ], 'mat_in', this_module, 'set_state -> state_vector')
         call zero_basis(state_out)
         do k = 1, kdim
            state_out(k)%state = mat_in(:,k)
         end do
      type is (state_matrix)
         call assert_shape(mat_in, [ N, N ], 'mat_in', this_module, 'set_state -> state_matrix')
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
      ! SVD
      real(wp)                           :: U_svd(rk,rk)
      real(wp)                           :: S_svd(rk)
      real(wp)                           :: V_svd(rk,rk)
      integer                            :: i, info
      character(len=128) :: msg

      if (size(U) < rk) then
         write(msg,'(A,I0)') 'Input krylov basis size incompatible with requested rank ', rk
         call stop_error(msg, module=this_module, procedure='generate_random_initial_condition')
      else
         call zero_basis(U)
         do i = 1,rk
            call U(i)%rand(.false.)
         end do
      end if
      call assert_shape(S, [ rk,rk ], 'S', this_module, 'generate_random_initial_condition')
      S = 0.0_wp

      ! perform QR
      allocate(Utmp(rk), source=U(:rk))
      call qr(Utmp, S, info)
      call check_info(info, 'qr', module=this_module, procedure='generate_random_initial_condition')
      ! perform SVD
      call svd(S(:rk,:rk), S_svd, U_svd, V_svd)
      S(:rk,:rk) = diag(S_svd)
      block
         class(abstract_vector_rdp), allocatable :: Xwrk(:)
         call linear_combination(Xwrk, Utmp, U_svd)
         call copy_basis(U, Xwrk)
      end block
      write(msg,'(A,I0,A,I0,A)') 'size(U) = [ ', size(U),' ]: filling the first ', rk, ' columns with noise.'
      call logger%log_information(msg, module=this_module, procedure='generate_random_initial_condition')
      return
   end subroutine

   !------------------------
   !-----     MISC     -----
   !------------------------

   function CALE(X,A,Q) result(res)
      real(wp), dimension(n,n) :: X, A, Q, res
      res = matmul(A, X) + matmul(X, transpose(A)) + Q
   end function CALE

   subroutine build_operator(A)
      !! Build the two-dimensional Laplace operator explicitly
      real(wp), intent(out) :: A(N,N)
      integer :: i, j, k

      A = -4.0_wp/dx2*eye(N)
      do i = 1, nx
         do j = 1, nx - 1
            k = (i-1)*nx + j
            A(k + 1, k) = 1.0_wp/dx2
            A(k, k + 1) = 1.0_wp/dx2
         end do 
      end do
      do i = 1, N-nx
         A(i, i + nx) = 1.0_wp/dx2
         A(i + nx, i) = 1.0_wp/dx2
      end do
      return
   end subroutine build_operator

   subroutine reconstruct_TQ(T, Q, A, D, E, tw)
      !! Reconstruct tridiagonal matrix T and orthogonal projector Q from dsytd2 output (A, D, E)
      real(wp), intent(out) :: T(N,N)
      real(wp), intent(out) :: Q(N,N)
      real(wp), intent(in)  :: A(N,N)
      real(wp), intent(in)  :: D(N)
      real(wp), intent(in)  :: E(N-1)
      real(wp), intent(in)  :: tw(N-1)

      ! internal variables
      real(wp)  :: Hi(N,N)
      real(wp)  :: vec(N,1)
      integer :: i

      ! Build orthogonal Q = H(1) @  H(2) @ ... @ H(n-1)
      Q = eye(N)
      do i = 1, N - 1
         vec          = 0.0_wp
         vec(i+1,1)   = 1.0_wp
         vec(i+2:N,1) = A(i+2:N,i)
         Hi           = eye(N) - tw(i) * matmul( vec, transpose(vec) )
         Q            = matmul( Q, Hi )
      end do

      ! Build tridiagonal T
      T = 0.0_wp
      do i = 1, N
         T(i,i) = D(i)
      end do
      do i = 1, N - 1
         T(i,i+1) = E(i)
         T(i+1,i) = E(i)
      end do

   end subroutine reconstruct_TQ

   subroutine pmat(A, name, row, col, p, prefix)
      real(wp), intent(in) :: A(:,:)
      character(len=*), optional, intent(in) :: name
      integer, optional, intent(in) :: row
      integer, optional, intent(in) :: col
      integer, optional, intent(in) :: p
      character(len=*), optional, intent(in) :: prefix
      ! internal
      integer :: i, row_, col_, p_, w
      character(len=128) :: name_, prefix_
      character(len=128) :: fmt

      row_ = optval(row, size(A,1))
      col_ = optval(col, size(A,2))
      p_   = optval(p, 5)
      prefix_ = optval(prefix, '')
      w = p_ + 3
      write(fmt,'("(A,3X,*(F",I0,".",I0,",1X))")') w, p_
      if (row_ .gt. size(A,1)) call logger%log_error('row > row(A)')
      if (col_ .gt. size(A,2)) call logger%log_error('col > col(A)')
      name_ = optval(name, 'mat')

      print *, name_
      do i = 1,row_
         print fmt, trim(prefix_), A(i,:col_)
      end do
      return
   end subroutine pmat

   subroutine pvec(U, name, col, p, prefix)
      class(abstract_vector_rdp), intent(in) :: U(:)
      character(len=*), optional, intent(in) :: name
      integer, optional, intent(in) :: col
      integer, optional, intent(in) :: p
      character(len=*), optional, intent(in) :: prefix
      ! internal
      real(wp), allocatable :: Umat(:,:)
      integer :: i, col_, kdim
      character(len=128) :: name_

      select type (U)
      type is (state_vector)
         kdim = size(U)
         col_ = optval(col, kdim)
         if (col_ .gt. kdim) call logger%log_error('col > col(U)')
         name_ = optval(name, 'mat')
         allocate(Umat(N,col_))
         call get_state(Umat, U(:col_), name)
         call pmat(Umat, name=name_, row=N, col=col_, p=p, prefix=prefix)
      end select
      return     
   end subroutine pvec

end module Laplacian2D_LTI_Lyapunov_Utils