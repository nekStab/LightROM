module LightROM_RiccatiUtils
   use LightKrylov
   use LightROM_AbstractLTIsystems
   use LightKrylov_utils, only : assert_shape
   implicit none

   ! scratch arrays
   real(kind=wp)         ,  allocatable   :: Swrk(:,:)

   private
   public :: apply_outerproduct_w, apply_p_outerproduct_w, precompute_NL

   !------------------------------
   !-----     INTERFACES     -----
   !------------------------------

   interface precompute_NL
      module procedure precompute_NL_K
      module procedure precompute_NL_S
   end interface

contains

   subroutine apply_outerproduct_w(Z, U, B, W)
      !! Computes the matrix product \( \mathbf{Z} = \mathbf{B} \mathbf{W} \mathbf{B}^T \mathbf{U} \) 
      class(abstract_vector),     intent(out)  :: Z(:)
      class(abstract_vector),     intent(in)   :: U(:)
      class(abstract_vector),     intent(in)   :: B(:)
      real(kind=wp),              intent(in)   :: W(:,:)
      ! internals
      integer                                  :: p, rk
      real(kind=wp),               allocatable :: wrk(:,:)

      p  = size(B)
      rk = size(U)
      allocate(wrk(1:p,1:rk)); wrk = 0.0_wp

      call assert_shape(W, (/ p, p /), 'apply_outerproduct_w', 'W')

      call mat_zero(Z)
      call mat_mult(wrk, B, U)
      call mat_mult(Z, B, matmul(W, wrk))
   
      return   
   end subroutine apply_outerproduct_w

   subroutine apply_p_outerproduct_w(M, UL, UR, B, W)
      !! Computes the matrix product \( \mathbf{M} = \mathbf{U}_L^T \mathbf{B} \mathbf{W} \mathbf{B}^T \mathbf{U}_R \) 
      real(kind=wp),              intent(out)  :: M(:,:)
      class(abstract_vector),     intent(in)   :: UL(:)
      class(abstract_vector),     intent(in)   :: UR(:)
      class(abstract_vector),     intent(in)   :: B(:)
      real(kind=wp),              intent(in)   :: W(:,:)
      ! internals
      real(kind=wp)                            :: BTUR(size(B),size(UR))
      real(kind=wp)                            :: ULTB(size(UL),size(B))

      call assert_shape(M, (/ size(UL), size(UR) /), 'apply_p_outerproduct_w', 'M')
      
      BTUR = 0.0_wp; ULTB = 0.0_wp; M = 0.0_wp
      call mat_mult(BTUR, B, UR)
      call mat_mult(ULTB, UL, B)
      
      M = matmul( ULTB, matmul( W, BTUR ) )
   
      return   
   end subroutine apply_p_outerproduct_w

   subroutine precompute_NL_K(N, X, K, B, W)
      !! Computes the matrix product \( \mathbf{N} = \mathbf{K} \mathbf{U}_L^T \mathbf{B} \mathbf{W} \mathbf{B}^T \mathbf{K} \) 
      class(abstract_vector),             intent(out)  :: N(:)
      class(abstract_sym_low_rank_state), intent(in)   :: X
      class(abstract_vector),             intent(in)   :: K(:)
      class(abstract_vector),             intent(in)   :: B(:)
      real(kind=wp),                      intent(in)   :: W(:,:)

      ! internals
      integer :: rk
      
      rk = size(K)
      if (.not.allocated(Swrk)) allocate(Swrk(1:rk,1:rk))
      Swrk = 0.0_wp

      call apply_p_outerproduct_w(Swrk, X%U, K, B, W)  ! (U.T) @ B @ R^(-1) @ B.T @ K
      call mat_mult(N, K, Swrk)                        ! K @ Swrk

      return
   end subroutine precompute_NL_K

   subroutine precompute_NL_S(N, X, U, B, W)
      !! Computes the matrix product \( \mathbf{N} = \mathbf{S} \mathbf{U}_L^T \mathbf{B} \mathbf{W} \mathbf{B}^T \mathbf{S} \) 
      real(kind=wp),                      intent(out)  :: N(:,:)
      class(abstract_sym_low_rank_state), intent(in)   :: X
      class(abstract_vector),             intent(in)   :: U(:)
      class(abstract_vector),             intent(in)   :: B(:)
      real(kind=wp),                      intent(in)   :: W(:,:)

      ! internals
      integer :: rk
      
      rk = size(U)
      if (.not.allocated(Swrk)) allocate(Swrk(1:rk,1:rk))
      Swrk = 0.0_wp

      call apply_p_outerproduct_w(Swrk, U, X%U, B, W)  !        U.T @ B @ R^(-1) @ B.T @ X%U
      N = matmul(X%S, matmul(Swrk, X%S))               ! X%S @ (U.T @ B @ R^(-1) @ B.T @ X%U) @ X%S

      return
   end subroutine precompute_NL_S

end module LightROM_RiccatiUtils