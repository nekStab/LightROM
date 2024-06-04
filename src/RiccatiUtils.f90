module LightROM_RiccatiUtils
   use LightKrylov
   use LightKrylov, only: dp
   use LightKrylov_AbstractVectors

   use LightROM_AbstractLTIsystems
   use LightKrylov_Utils, only : assert_shape
   implicit none

   ! scratch arrays
   real(dp),       allocatable   :: Swrk(:,:)

   private

   public :: apply_outerprod_w
   public :: apply_premult_outerprod_w
   public :: precompute_NL

   !------------------------------
   !-----     INTERFACES     -----
   !------------------------------

   interface apply_outerprod_w
      module procedure apply_outerprod_matrix_w_rdp
   end interface

   interface apply_premult_outerprod_w
      module procedure apply_premult_outerprod_matrix_w_rdp
   end interface

   interface precompute_NL
      module procedure precompute_NL_K_rdp
      module procedure precompute_NL_S_rdp
   end interface

contains

   subroutine apply_outerprod_matrix_w_rdp(Z, U, B, W)
      !! Computes the matrix product \( \mathbf{Z} = \mathbf{B} \mathbf{W} \mathbf{B}^T \mathbf{U} \) 
      class(abstract_vector_rdp), intent(out)  :: Z(:)
      class(abstract_vector_rdp), intent(in)   :: U(:)
      class(abstract_vector_rdp), intent(in)   :: B(:)
      real(dp),                   intent(in)   :: W(:,:)
      ! internals
      integer                                  :: p, rk, i
      real(dp),                    allocatable :: wrk(:,:)

      p  = size(B)
      rk = size(U)
      allocate(wrk(1:p,1:rk)); wrk = 0.0_dp

      call assert_shape(W, (/ p, p /), 'apply_outerprod_matrix_w_rdp', 'W')

      do i = 1, size(Z)
         call Z(i)%zero()
      end do
      call innerprod_matrix(wrk, B, U)
      block
         class(abstract_vector_rdp), allocatable :: Xwrk(:)
         call linear_combination(Xwrk, B, matmul(W, wrk))
         call copy_basis(Z, Xwrk)
      end block
   
      return   
   end subroutine apply_outerprod_matrix_w_rdp

   subroutine apply_premult_outerprod_matrix_w_rdp(M, UL, UR, B, W)
      !! Computes the matrix product \( \mathbf{M} = \mathbf{U}_L^T \mathbf{B} \mathbf{W} \mathbf{B}^T \mathbf{U}_R \) 
      real(dp),                   intent(out)  :: M(:,:)
      class(abstract_vector_rdp), intent(in)   :: UL(:)
      class(abstract_vector_rdp), intent(in)   :: UR(:)
      class(abstract_vector_rdp), intent(in)   :: B(:)
      real(dp),                   intent(in)   :: W(:,:)
      ! internals
      real(dp)                                 :: BTUR(size(B),size(UR))
      real(dp)                                 :: ULTB(size(UL),size(B))

      call assert_shape(M, (/ size(UL), size(UR) /), 'apply_premult_outerprod_matrix_w_rdp', 'M')
      
      BTUR = 0.0_dp; ULTB = 0.0_dp; M = 0.0_dp
      call innerprod_matrix(BTUR, B, UR)
      call innerprod_matrix(ULTB, UL, B)
      
      M = matmul( ULTB, matmul( W, BTUR ) )
   
      return   
   end subroutine apply_premult_outerprod_matrix_w_rdp

   subroutine precompute_NL_K_rdp(N, X, K, B, W)
      !! Computes the matrix product \( \mathbf{N} = \mathbf{K} \mathbf{U}_L^T \mathbf{B} \mathbf{W} \mathbf{B}^T \mathbf{K} \) 
      class(abstract_vector_rdp),             intent(out)  :: N(:)
      class(abstract_sym_low_rank_state_rdp), intent(in)   :: X
      class(abstract_vector_rdp),             intent(in)   :: K(:)
      class(abstract_vector_rdp),             intent(in)   :: B(:)
      real(dp),                               intent(in)   :: W(:,:)

      ! internals
      integer :: rk
      
      rk = size(K)
      if (.not.allocated(Swrk)) allocate(Swrk(1:rk,1:rk))
      Swrk = 0.0_dp

      call apply_premult_outerprod_w(Swrk, X%U, K, B, W)  ! (U.T) @ B @ R^(-1) @ B.T @ K
      block
         class(abstract_vector_rdp), allocatable :: Xwrk(:)
         call linear_combination(Xwrk, K, Swrk)                 ! K @ Swrk
         call copy_basis(N, Xwrk)
      end block

      return
   end subroutine precompute_NL_K_rdp

   subroutine precompute_NL_S_rdp(N, X, U, B, W)
      !! Computes the matrix product \( \mathbf{N} = \mathbf{S} \mathbf{U}_L^T \mathbf{B} \mathbf{W} \mathbf{B}^T \mathbf{S} \) 
      real(dp),                               intent(out)  :: N(:,:)
      class(abstract_sym_low_rank_state_rdp), intent(in)   :: X
      class(abstract_vector_rdp),             intent(in)   :: U(:)
      class(abstract_vector_rdp),             intent(in)   :: B(:)
      real(dp),                               intent(in)   :: W(:,:)

      ! internals
      integer :: rk
      
      rk = size(U)
      if (.not.allocated(Swrk)) allocate(Swrk(1:rk,1:rk))
      Swrk = 0.0_dp

      call apply_premult_outerprod_w(Swrk, U, X%U, B, W)  !        U.T @ B @ R^(-1) @ B.T @ X%U
      N = matmul(X%S, matmul(Swrk, X%S))                                   ! X%S @ (U.T @ B @ R^(-1) @ B.T @ X%U) @ X%S

      return
   end subroutine precompute_NL_S_rdp

end module LightROM_RiccatiUtils