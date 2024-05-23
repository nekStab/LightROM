module Laplacian2D_LTI_Riccati_Operators
   use Laplacian2D_LTI_Riccati_Base
   !> LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov_utils
   !> Standard Library.
   use stdlib_math, only : linspace
   use stdlib_optval, only : optval
   use stdlib_linalg, only : eye
   implicit none

   ! operator
   public :: CARE, build_operator, laplacian, laplacian_mat

   !-----------------------------------
   !-----     LAPLACE OPERATOR    -----
   !-----------------------------------

   type, extends(abstract_linop), public :: laplace_operator
   contains
      private
      procedure, pass(self), public :: matvec  => direct_matvec_laplace
      procedure, pass(self), public :: rmatvec => direct_matvec_laplace     ! dummy since Lyapunov equation for Laplacian is symmetric
   end type laplace_operator

contains

   function CARE(X,A,Q,BRinvBT) result(Y)
      real(kind=wp), dimension(n,n) :: X, A, Q, BRinvBT, Y
      Y = matmul(transpose(A), X) + matmul(X, A) + Q - matmul(X, matmul(BRinvBT, X))
   end function CARE

   !-----     TYPE-BOUND PROCEDURE FOR LAPLACE OPERATOR    -----

   subroutine direct_matvec_laplace(self, vec_in, vec_out)
      !> Linear Operator.
      class(laplace_operator),intent(in)  :: self
      !> Input vector.
      class(abstract_vector) , intent(in)  :: vec_in
      !> Output vector.
      class(abstract_vector) , intent(out) :: vec_out
      select type(vec_in)
      type is (state_vector)
         select type(vec_out)
         type is (state_vector)
            call laplacian(vec_out%state, vec_in%state)
         end select
      end select
      return
   end subroutine direct_matvec_laplace

   !---------------------------
   !-----    Laplacian    -----
   !---------------------------

   subroutine build_operator(A)
      !! Build the two-dimensional Laplace operator explicitly
      real(kind=wp), intent(out) :: A(N,N)
      integer i, j, k

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

   subroutine laplacian(vec_out, vec_in)
      
      !> State vector.
      real(kind=wp)  , dimension(:), intent(in)  :: vec_in
      !> Time-derivative.
      real(kind=wp)  , dimension(:), intent(out) :: vec_out

      !> Internal variables.
      integer             :: i, j, in
      
      in = 1
      vec_out(in)       = (                                  - 4*vec_in(in) + vec_in(in + 1) + vec_in(in + nx)) / dx2
      do in = 2, nx - 1
         vec_out(in)    = (                   vec_in(in - 1) - 4*vec_in(in) + vec_in(in + 1) + vec_in(in + nx)) / dx2
      end do
      in = nx
      vec_out(in)       = (                   vec_in(in - 1) - 4*vec_in(in)                  + vec_in(in + nx)) / dx2
      !
      do i = 2, nx-1
         in = (i-1)*nx + 1
         vec_out(in)    = ( vec_in(in - nx)                  - 4*vec_in(in) + vec_in(in + 1) + vec_in(in + nx)) / dx2
         do j = 2, nx - 1
            in = (i-1)*nx + j
            vec_out(in) = ( vec_in(in - nx) + vec_in(in - 1) - 4*vec_in(in) + vec_in(in + 1) + vec_in(in + nx)) / dx2 
         end do
         in = (i-1)*nx + nx
         vec_out(in)    = ( vec_in(in - nx) + vec_in(in - 1) - 4*vec_in(in)                  + vec_in(in + nx)) / dx2
      end do
      !
      in = N - nx + 1
      vec_out(in)       = ( vec_in(in - nx)                  - 4*vec_in(in) + vec_in(in + 1)                  ) / dx2
      do in = N - nx + 2, N - 1
         vec_out(in)    = ( vec_in(in - nx) + vec_in(in - 1) - 4*vec_in(in) + vec_in(in + 1)                  ) / dx2
      end do
      in = N
      vec_out(in)       = ( vec_in(in - nx) + vec_in(in - 1) - 4*vec_in(in)                                   ) / dx2
         
      return
   end subroutine laplacian

   subroutine laplacian_mat(flat_mat_out, flat_mat_in, transpose)
   
      !> State vector.
      real(kind=wp)  , dimension(:), intent(in)  :: flat_mat_in
      !> Time-derivative.
      real(kind=wp)  , dimension(:), intent(out) :: flat_mat_out
      !> Transpose
      logical, optional :: transpose
      logical           :: trans
      
      !> Internal variables.
      integer :: j
      real(kind=wp), dimension(N,N) :: mat, dmat
      
      !> Deal with optional argument
      trans = optval(transpose,.false.)
      
      !> Sets the internal variables.
      mat  = reshape(flat_mat_in(1:N**2),(/N, N/))
      dmat = 0.0_wp
      
      if (trans) then
          do j = 1,N
             call laplacian(dmat(j,:), mat(j,:))
          end do
      else
          do j = 1,N
             call laplacian(dmat(:,j), mat(:,j))
          end do
      endif

      !> Reshape for output
      flat_mat_out = reshape(dmat, shape(flat_mat_in))
       
      return
   end subroutine laplacian_mat

end module Laplacian2D_LTI_Riccati_Operators