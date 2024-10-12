module Laplacian2D_LTI_Riccati_Base
   ! Standard Library.
   use stdlib_optval, only : optval
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov, only : wp => dp
   use LightKrylov_Logger
   use LightKrylov_Utils, only : assert_shape
   use LightKrylov_AbstractVectors
   ! LightROM
   use LightROM_AbstractLTIsystems ! LR_state
   implicit none

   private :: this_module
   character*128, parameter :: this_module = 'Laplacian2D_LTI_Riccati_Base'
   ! problem parameters
   public  :: N, nx, dx, dx2, L, rk_b, rk_c
   ! problem definition
   public  :: B, CT, Qc, Rinv
   ! derived data
   public  :: Bdata, CTdata, CTQcC, BRinvBTdata, CTQcCdata

   !------------------------------
   !-----     PARAMETERS     -----
   !------------------------------

   ! --> Mesh related parameters.
   real(wp), parameter :: L  = 1.0_wp  !> Domain length
   integer,  parameter :: nx = 4      !> Number of grid points per direction
   integer,  parameter :: N  = nx**2   !> total number of grid points
   real(wp), parameter :: dx = L/nx    !> Grid size.
   real(wp), parameter :: dx2= dx**2   !> Grid size.
   integer,  parameter :: rk_b = 1     !> rank of the RHS
   integer,  parameter :: rk_c = 1     !> rank of Q = CTC

   !-------------------------------------------------------
   !-----     LIGHTKRYLOV SYM LOW RANK STATE TYPE     -----
   !-------------------------------------------------------

   type, extends(abstract_sym_low_rank_state_rdp), public :: LR_state
   contains
      private
      procedure, pass(self), public :: initialize_LR_state
   end type LR_state


   !-------------------------------------------
   !-----     LIGHTKRYLOV VECTOR TYPE     -----
   !-------------------------------------------

   type, extends(abstract_vector_rdp), public :: state_vector
      real(wp) :: state(N) = 0.0_wp
      contains
      private
      procedure, pass(self), public :: zero => vector_zero
      procedure, pass(self), public :: dot => vector_dot
      procedure, pass(self), public :: scal => vector_scal
      procedure, pass(self), public :: axpby => vector_axpby
      procedure, pass(self), public :: rand => vector_rand
      procedure, pass(self), public :: get_size => vector_get_size
   end type state_vector

   !-------------------------------------------
   !-----     LIGHTKRYLOV VECTOR TYPE     -----
   !-------------------------------------------

   type, extends(abstract_vector_rdp), public :: state_matrix
      real(wp) :: state(N**2) = 0.0_wp
   contains
      private
      procedure, pass(self), public :: zero => matrix_zero
      procedure, pass(self), public :: dot => matrix_dot
      procedure, pass(self), public :: scal => matrix_scal
      procedure, pass(self), public :: axpby => matrix_axpby
      procedure, pass(self), public :: rand => matrix_rand
      procedure, pass(self), public :: get_size => matrix_get_size
   end type state_matrix

   type(state_vector)       :: B(rk_b)
   type(state_vector)       :: CT(rk_c)
   real(wp)                 :: Qc(rk_c,rk_c)
   real(wp)                 :: Rinv(rk_b,rk_b)

   real(wp)                 :: Bdata(N,rk_b)
   real(wp)                 :: CTdata(N,rk_c)
   real(wp)                 :: CTQcC(N**2)
   real(wp)                 :: CTQcCdata(N,N)
   real(wp)                 :: BRinvBTdata(N,N)

contains

   !-----     TYPE-BOUND PROCEDURE FOR VECTORS     -----

   subroutine vector_zero(self)
      class(state_vector), intent(inout) :: self
      self%state = 0.0_wp
      return
   end subroutine vector_zero

   real(wp) function vector_dot(self, vec) result(alpha)
      class(state_vector),        intent(in) :: self
      class(abstract_vector_rdp), intent(in) :: vec
      select type(vec)
      type is (state_vector)
         alpha = dot_product(self%state, vec%state)
      end select
      return
   end function vector_dot

   integer function vector_get_size(self) result(N)
     class(state_vector), intent(in) :: self
     N = nx
     return
   end function vector_get_size

   subroutine vector_scal(self, alpha)
      class(state_vector), intent(inout) :: self
      real(wp),            intent(in)    :: alpha
      self%state = self%state * alpha
      return
   end subroutine vector_scal

   subroutine vector_axpby(self, alpha, vec, beta)
      class(state_vector),        intent(inout) :: self
      class(abstract_vector_rdp), intent(in)    :: vec
      real(wp),                   intent(in)    :: alpha, beta
      select type(vec)
      type is (state_vector)
         self%state = alpha*self%state + beta*vec%state
      end select
      return
   end subroutine vector_axpby

   subroutine vector_rand(self, ifnorm)
      class(state_vector), intent(inout) :: self
      logical, optional,   intent(in)    :: ifnorm
      ! internals
      logical :: normalize
      real(wp) :: alpha
      normalize = optval(ifnorm,.true.)
      call random_number(self%state)
      if (normalize) then
         alpha = self%norm()
         call self%scal(1.0/alpha)
      endif
      return
   end subroutine vector_rand

   !-----     TYPE-BOUND PROCEDURE FOR MATRICES     -----

   subroutine matrix_zero(self)
      class(state_matrix), intent(inout) :: self
      self%state = 0.0_wp
      return
   end subroutine matrix_zero

   real(wp) function matrix_dot(self, vec) result(alpha)
      class(state_matrix),        intent(in) :: self
      class(abstract_vector_rdp), intent(in) :: vec
      select type(vec)
      type is(state_matrix)
          alpha = dot_product(self%state, vec%state)
      end select
      return
   end function matrix_dot

   integer function matrix_get_size(self) result(N)
     class(state_matrix), intent(in) :: self
     N = N
     return
   end function matrix_get_size

   subroutine matrix_scal(self, alpha)
      class(state_matrix), intent(inout) :: self
      real(wp),            intent(in)    :: alpha
      self%state = self%state * alpha
      return
   end subroutine matrix_scal  

   subroutine matrix_axpby(self, alpha, vec, beta)
      class(state_matrix),        intent(inout) :: self
      class(abstract_vector_rdp), intent(in)    :: vec
      real(wp),                   intent(in)    :: alpha, beta
      select type(vec)
      type is(state_matrix)
          self%state = alpha*self%state + beta*vec%state
      end select
      return
   end subroutine matrix_axpby

   subroutine matrix_rand(self, ifnorm)
      class(state_matrix), intent(inout) :: self
      logical, optional,   intent(in)    :: ifnorm
      ! internals
      logical :: normalize
      real(wp)      :: alpha
      normalize = optval(ifnorm, .true.)
      call random_number(self%state)
      if (normalize) then
         alpha = self%norm()
         call self%scal(1.0/alpha)
      endif
      return
   end subroutine matrix_rand  

   !-----------------------------------------------------------------------
   !-----     TYPE BOUND PROCEDURE FOR SYM LOW RANK REPRESENTATION    -----
   !-----------------------------------------------------------------------

   subroutine initialize_LR_state(self, U, S, rk, rkmax)
      class(LR_state),            intent(inout) :: self
      class(abstract_vector_rdp), intent(in)    :: U(:)
      real(wp),                   intent(in)    :: S(:,:)
      integer,                    intent(in)    :: rk
      integer, optional,          intent(in)    :: rkmax

      ! internals
      real(wp), allocatable :: R(:, :)
      integer :: i, n, rka, info

      n = size(U)
      call assert_shape(S, [n,n], 'S', this_module, 'initialize_LR_state')

      ! optional size argument
      if (present(rkmax)) then
         self%rk = rkmax - 1
         rka = rkmax
      else
         self%rk = rk
         rka = rk + 1
      end if

      select type (U)
      type is (state_vector)
         ! allocate & initialize
         allocate(self%U(rka), source=U(1)); call zero_basis(self%U)
         allocate(self%S(rka,rka)); self%S = 0.0_wp
         ! copy inputs
         if (self%rk > n) then   ! copy the full IC into self%U
            call copy_basis(self%U(1:n), U)
            self%S(1:n,1:n) = S
         else  ! fill the first self%rk columns of self%U with the first self%rk columns of the IC
            call copy_basis(self%U(1:self%rk), U(1:self%rk))
            self%S(1:self%rk,1:self%rk) = S(1:self%rk,1:self%rk)
         end if
         ! top up basis (to rka for rank-adaptivity) with orthonormal columns if needed
         if (rka > n) then
            do i = n+1, rka
               call self%U(i)%rand()
            end do
            allocate(R(rka,rka)); R = 0.0_wp
            call qr(self%U, R, info)
            call check_info(info, 'qr', module=this_module, procedure='initialize_LR_state')
         end if
      end select
      return
   end subroutine initialize_LR_state

end module Laplacian2D_LTI_Riccati_Base