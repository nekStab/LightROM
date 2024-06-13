module LightROM_RiccatiUtils
   use LightKrylov
   use LightKrylov, only: wp => dp
   use LightKrylov_AbstractVectors

   use LightROM_AbstractLTIsystems
   use LightKrylov_Utils, only : assert_shape
   implicit none

   private

   public :: apply_outerprod_w
   public :: apply_premult_outerprod_w
   public :: precompute_NL

   !------------------------------
   !-----     INTERFACES     -----
   !------------------------------

   interface apply_outerprod_w
      module procedure apply_outerprod_w_vector_rdp
      module procedure apply_outerprod_w_basis_rdp
   end interface

   interface apply_premult_outerprod_w
      module procedure apply_premult_outerprod_w_vector_rdp
      module procedure apply_premult_outerprod_w_basis_rdp
   end interface

   interface precompute_NL
      module procedure precompute_NL_K_rdp
      module procedure precompute_NL_S_rdp
   end interface

contains

   subroutine apply_outerprod_w_vector_rdp(z, u, B, W)
      !! Computes the matrix product \( \mathbf{z} = \mathbf{B} \mathbf{W} \mathbf{B}^T \mathbf{u} \) 
      class(abstract_vector_rdp), intent(out)  :: z
      class(abstract_vector_rdp), intent(in)   :: u
      class(abstract_vector_rdp), i      !call linear_combination(Utmp, X%U, Swrk(1:rk,1:rk))ntent(in)   :: B(:)
      real(wp),                   intent(in)   :: W(:,:)
      ! internals
      real(wp)                                 :: wrk(size(B))

      call assert_shape(W, (/ size(B), size(B) /), 'apply_outerprod_w_vector_rdp', 'W')

      call innerprod(wrk, B, u)
      block
         class(abstract_vector_rdp), allocatable :: xwrk
         call linear_combination(xwrk, B, matmul(W, wrk))
         call z%zero(); call z%add(xwrk)
      end block

      return   
   end subroutine apply_outerprod_w_vector_rdp

   subroutine apply_outerprod_w_basis_rdp(Z, U, B, W)
      !! Computes the matrix product \( \mathbf{Z} = \mathbf{B} \mathbf{W} \mathbf{B}^T \mathbf{U} \) 
      class(abstract_vector_rdp), intent(out)  :: Z(:)
      class(abstract_vector_rdp), intent(in)   :: U(:)
      class(abstract_vector_rdp), intent(in)   :: B(:)
      real(wp),                   intent(in)   :: W(:,:)
      ! internals
      real(wp)                                 :: wrk(size(B),size(U))

      call assert_shape(W, (/ size(B), size(B) /), 'apply_outerprod_w_basis_rdp', 'W')

      call zero_basis(Z)
      call innerprod(wrk, B, U)
      block
         class(abstract_vector_rdp), allocatable :: Xwrk(:)
         call linear_combination(Xwrk, B, matmul(W, wrk))
         call copy_basis(Z, Xwrk)
      end block
   
      return   
   end subroutine apply_outerprod_w_basis_rdp

   subroutine apply_premult_outerprod_w_vector_rdp(m, uL, uR, B, W)
      !! Computes the matrix product \( \mathbf{M} = \mathbf{U}_L^T \mathbf{B} \mathbf{W} \mathbf{B}^T \mathbf{U}_R \) 
      real(wp),                   intent(out)  :: m
      class(abstract_vector_rdp), intent(in)   :: uL
      class(abstract_vector_rdp), intent(in)   :: uR
      class(abstract_vector_rdp), intent(in)   :: B(:)
      real(wp),                   intent(in)   :: W(:,:)
      ! internals
      real(wp)                                 :: BTuR(size(B))
      real(wp)                                 :: uLTB(size(B))
      
      call assert_shape(W, (/ size(B), size(B) /), 'apply_premult_outerprod_w_vector_rdp', 'W')

      BTuR = 0.0_wp; uLTB = 0.0_wp; m = 0.0_wp
      call innerprod(BTuR, B, uR)
      call innerprod(uLTB, B, uL)
      
      m = dot_product( uLTB, matmul( W, BTuR ) )
   
      return   
   end subroutine apply_premult_outerprod_w_vector_rdp

   subroutine apply_premult_outerprod_w_basis_rdp(M, UL, UR, B, W)
      !! Computes the matrix product \( \mathbf{M} = \mathbf{U}_L^T \mathbf{B} \mathbf{W} \mathbf{B}^T \mathbf{U}_R \) 
      real(wp),                   intent(out)  :: M(:,:)
      class(abstract_vector_rdp), intent(in)   :: UL(:)
      class(abstract_vector_rdp), intent(in)   :: UR(:)
      class(abstract_vector_rdp), intent(in)   :: B(:)
      real(wp),                   intent(in)   :: W(:,:)
      ! internals
      real(wp)                                 :: BTUR(size(B),size(UR))
      real(wp)                                 :: ULTB(size(UL),size(B))

      call assert_shape(W, (/ size(B),  size(B)  /), 'apply_premult_outerprod_w_basis_rdp', 'W')
      call assert_shape(M, (/ size(UL), size(UR) /), 'apply_premult_outerprod_w_basis_rdp', 'M')
      
      BTUR = 0.0_wp; ULTB = 0.0_wp; M = 0.0_wp
      call innerprod(BTUR, B, UR)
      call innerprod(ULTB, UL, B)
      
      M = matmul( ULTB, matmul( W, BTUR ) )
   
      return   
   end subroutine apply_premult_outerprod_w_basis_rdp

   subroutine precompute_NL_K_rdp(N, X, K, B, W)
      !! Computes the matrix product \( \mathbf{N} = \mathbf{K} \mathbf{U}_L^T \mathbf{B} \mathbf{W} \mathbf{B}^T \mathbf{K} \) 
      class(abstract_vector_rdp),             intent(out)  :: N(:)
      class(abstract_sym_low_rank_state_rdp), intent(in)   :: X
      class(abstract_vector_rdp),             intent(in)   :: K(:)
      class(abstract_vector_rdp),             intent(in)   :: B(:)
      real(wp),                               intent(in)   :: W(:,:)

      ! internals
      real(wp)                                             :: wrk(X%rk,X%rk)

      call assert_shape(W, (/ size(B), size(B) /), 'precompute_NL_K_rdp', 'W')
      
      call apply_premult_outerprod_w(wrk, X%U(1:X%rk), K, B, W) ! (U.T) @ B @ R^(-1) @ B.T @ K
      block
         class(abstract_vector_rdp), allocatable :: Xwrk(:)
         call linear_combination(Xwrk, K, wrk)                  ! K @ (U.T @ B @ R^(-1) @ B.T @ K)
         call copy_basis(N, Xwrk)
      end block

      return
   end subroutine precompute_NL_K_rdp

   subroutine precompute_NL_S_rdp(N, X, U, B, W)
      !! Computes the matrix product \( \mathbf{N} = \mathbf{S} \mathbf{U}_L^T \mathbf{B} \mathbf{W} \mathbf{B}^T \mathbf{S} \) 
      real(wp),                               intent(out)  :: N(:,:)
      class(abstract_sym_low_rank_state_rdp), intent(in)   :: X
      class(abstract_vector_rdp),             intent(in)   :: U(:)
      class(abstract_vector_rdp),             intent(in)   :: B(:)
      real(wp),                               intent(in)   :: W(:,:)

      ! internals
      real(wp)                                             :: wrk(X%rk,X%rk)

      call assert_shape(W, (/ size(B), size(B) /), 'precompute_NL_S_rdp', 'W')

      call apply_premult_outerprod_w(wrk, U, X%U(1:X%rk), B, W)       !        U.T @ B @ R^(-1) @ B.T @ X%U
      N = matmul(X%S(1:X%rk,1:X%rk), matmul(wrk, X%S(1:X%rk,1:X%rk))) ! X%S @ (U.T @ B @ R^(-1) @ B.T @ X%U) @ X%S

      return
   end subroutine precompute_NL_S_rdp

end module LightROM_RiccatiUtils