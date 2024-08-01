module Laplacian2D_LTI_Lyapunov_Utils
   ! Standard Library.
   use stdlib_math, only : linspace
   use stdlib_optval, only : optval
   use stdlib_linalg, only : eye, diag, svd, svdvals
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

   character*128, parameter :: this_module = 'Laplacian2D_LTI_Lyapunov_Utils'

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
      call svd(S(:,1:rk), S_svd(1:rk), U_svd(:,1:rk), V_svd(1:rk,1:rk))
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

   function CALE(X, A, Q, iftrans) result(Y)
      real(wp), dimension(n,n) :: X, A, Q, Y
      logical, optional :: iftrans
      logical           :: trans

      trans = optval(iftrans, .false.)

      if (trans) then
         Y = matmul(transpose(A), X) + matmul(X, A) + Q
      else
         Y = matmul(A, X) + matmul(X, transpose(A)) + Q
      end if
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

   subroutine print_svdvals(X_out, tag)
      !! Compute singular values and print the non-zero ones
      real(wp),         intent(in) :: X_out(:,:)
      character(len=*), intent(in) :: tag
      ! internal
      real(wp), allocatable :: svals(:)
      integer i, n, rki

      N = size(X_out)
      svals = svdvals(X_out)
      rki = N
      do i = 1, N
         if (svals(i) < 1e-14) then
             rki = i
             exit
         end if
      end do
      write(*,'(A,1X,A4,1X,I3,1X)', ADVANCE='NO') '  OUTPUT_SVD', trim(tag), rki
      do i = 1, rki
         write(*,'(E8.2,1x)', ADVANCE='NO') svals(i)
      end do
      write(*,*) 
      return
   end subroutine print_svdvals

   subroutine print_info(if_rkad, if_kexpm, Tend, dt, rk, TO)
      logical,  intent(in) :: if_rkad
      logical,  intent(in) :: if_kexpm
      real(wp), intent(in) :: Tend
      real(wp), intent(in) :: dt
      integer,  intent(in) :: rk
      integer,  intent(in) :: TO
      write(*,*) 
      write(*,*) "!----------------------------"
      write(*,*) 
      write(*,*) "DLRA Lyapunov equation:"
      write(*,'(A,F12.6)') "    Tend  = ", Tend
      write(*,'(A,E12.6)') "    dt    = ", dt
      write(*,*) "   rk0   = ", rk
      write(*,*) "   order = ", TO
      write(*,*) 
      write(*,*) "Parameters:"
      write(*,*) "   rank_adaptive = ", if_rkad
      write(*,*) "   Krylov-exptA  = ", if_kexpm
      write(*,*) 
      write(*,*) "!----------------------------"
      write(*,*) 
      return
   end subroutine print_info

end module Laplacian2D_LTI_Lyapunov_Utils
