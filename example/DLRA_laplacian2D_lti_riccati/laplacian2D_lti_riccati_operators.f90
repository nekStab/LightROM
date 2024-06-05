module Laplacian2D_LTI_Riccati_Operators
   ! Standard Library.
   use stdlib_optval, only : optval
   use stdlib_linalg, only : eye
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov, only : wp => dp
   use LightKrylov_Utils ! svd
   ! LightROM
   use LightROM_AbstractLTIsystems ! abstract_lti_system
   ! Laplacian
   use Laplacian2D_LTI_Riccati_Base
   
   implicit none

   ! operator
   public :: build_operator, laplacian, laplacian_mat, exptA

   !-----------------------------------------------
   !-----     LIGHTKRYLOV LTI SYSTEM TYPE     -----
   !-----------------------------------------------

   type, extends(abstract_lti_system_rdp), public :: lti_system
   contains
      private
      procedure, pass(self), public :: initialize_lti_system
   end type lti_system

   !-----------------------------------
   !-----     LAPLACE OPERATOR    -----
   !-----------------------------------

   type, extends(abstract_linop_rdp), public :: laplace_operator
   contains
      private
      procedure, pass(self), public :: matvec  => direct_matvec_laplace
      procedure, pass(self), public :: rmatvec => direct_matvec_laplace     ! dummy since Lyapunov equation for Laplacian is symmetric
   end type laplace_operator

contains

   !-----     TYPE-BOUND PROCEDURE FOR LAPLACE OPERATOR    -----

   subroutine direct_matvec_laplace(self, vec_in, vec_out)
      !> Linear Operator.
      class(laplace_operator),    intent(in)  :: self
      !> Input vector.
      class(abstract_vector_rdp), intent(in)  :: vec_in
      !> Output vector.
      class(abstract_vector_rdp), intent(out) :: vec_out
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
      real(wp), intent(out) :: A(N,N)
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
      real(wp), dimension(:), intent(in)  :: vec_in
      !> Time-derivative.
      real(wp), dimension(:), intent(out) :: vec_out

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
      real(wp), dimension(:), intent(in)  :: flat_mat_in
      !> Time-derivative.
      real(wp), dimension(:), intent(out) :: flat_mat_out
      !> Transpose
      logical, optional :: transpose
      logical           :: trans
      
      !> Internal variables.
      integer :: j
      real(wp), dimension(N,N) :: mat, dmat
      
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

   !--------------------------------------
   !-----     EXP(tA) SUBROUTINE     -----
   !--------------------------------------

   subroutine exptA(vec_out, A, vec_in, tau, info, trans)
      !! Subroutine for the exponential propagator that conforms with the abstract interface
      !! defined in expmlib.f90
      class(abstract_vector_rdp),  intent(out)   :: vec_out
      !! Output vector
      class(abstract_linop_rdp),   intent(inout) :: A
      !! Linear operator
      class(abstract_vector_rdp),  intent(in)    :: vec_in
      !! Input vector.
      real(wp),                    intent(in)    :: tau
      !! Integration horizon
      integer,                     intent(out)   :: info
      !! Information flag
      logical, optional,           intent(in)    :: trans
      logical                                    :: transpose
      !! Direct or Adjoint?

      ! optional argument
      transpose = optval(trans, .false.)

      ! time integrator
      select type (vec_in)
      type is (state_vector)
         select type (vec_out)
         type is (state_vector)
            select type (A)
            type is (laplace_operator)
               call k_exptA(vec_out, A, vec_in, tau, info, transpose)
            end select
         end select
      end select

   end subroutine exptA

   !--------------------------------------------------------
   !-----     TYPE BOUND PROCEDURES FOR LTI SYSTEMS    -----
   !--------------------------------------------------------

   subroutine initialize_lti_system(self, A, B, CT, D)
      class(lti_system),           intent(inout) :: self
      class(abstract_linop_rdp),   intent(in)    :: A
      class(abstract_vector_rdp),  intent(in)    :: B(:)
      class(abstract_vector_rdp),  intent(in)    :: CT(:)
      real(wp),          optional, intent(in)    :: D(:,:)

      ! internal variables
      integer                                :: rk_b, rk_c

      ! Operator
      select type (A)
      type is (laplace_operator)
         allocate(self%A, source=A)
      end select
      ! Input
      select type (B)
      type is (state_vector)
         rk_b = size(B)
         allocate(self%B(1:rk_b), source=B(1:rk_b))
      end select
      ! Output
      select type (CT)
         type is (state_vector)
         rk_c = size(CT)
         allocate(self%CT(1:rk_c), source=CT(1:rk_c))
      end select
      ! Throughput
      allocate(self%D(1:rk_c, 1:rk_b))
      if (present(D)) then
         call assert_shape(D, (/ rk_c, rk_b /), 'initialize_lti_system', 'D')
         self%D = D
      else
         self%D = 0.0_wp
      end if
      return
   end subroutine initialize_lti_system

end module Laplacian2D_LTI_Riccati_Operators