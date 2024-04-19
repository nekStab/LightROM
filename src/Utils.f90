module LightROM_utils
   use LightKrylov
   use LightKrylov_utils
   use LightROM_AbstractLTIsystems

   use stdlib_linalg, only : eye, diag
   use stdlib_optval, only : optval
   implicit none 

   private
   public Approximate_Balanced_Truncation, ROM_Petrov_Galerkin_Projection, ROM_Galerkin_Projection

contains

   subroutine approximate_Balanced_Truncation(T,S,ST,X,Y)
      !! Computes the the biorthogonal balancing transformation \( \mathbf{T}, \mathbf{S}^T \) from the
      !! low-rank approximations of the obervability and controlability Gramians, \( \mathbf{W}_o \) and 
      !! \( \mathbf{W}_c \) respectively, given as:
      !! \[ 
      !!    \mathbf{W}_o = \mathbf{X}_o \mathbf{S}_o \mathbf{X}_o^T 
      !!    \quad \text{and} \quad 
      !!    \mathbf{W}_c = \mathbf{Y}_c \mathbf{S}_c \mathbf{Y}_c^T
      !! \]
      !!
      !! Given the SVD of the cross-Gramians:
      !! $$ \mathbf{S}_c^T \mathbf{Y}_c^T \mathbf{X}_o \mathbf{S}_o = \mathbf{U} \mathbf{S} \mathbf{V}^T $$
      !! the balancing transformation and its inverse are given by:
      !! \[ \begin{align}
      !!            \mathbf{T}   &= \mathbf{X}_o \mathbf{S}_o^{1/2} \mathbf{V} \mathbf{S}^{-1/2} \\
      !!            \mathbf{S}^T &= \mathbf{Y}_c \mathbf{S}_c^{1/2} \mathbf{U} \mathbf{S}^{-1/2} 
      !! \end{align} \]
      !! Note: In the current implementation, the numerical rank of the SVD is not considered.
      class(abstract_vector),              intent(out)   :: T(:)
      !! Balancing transformation
      real(kind=wp),                       intent(out)   :: S(:)
      !! Singular values of the BT
      class(abstract_vector),              intent(out)   :: ST(:)
      !! Inverse balancing transformation
      class(abstract_sym_low_rank_state),  intent(inout) :: X
      !! Low-rank representation of the Controllability Gramian
      class(abstract_sym_low_rank_state),  intent(inout) :: Y
      !! Low-rank representation of the Observability Gramian

      ! internal variables
      integer                                :: i, rkc, rko, rk, rkmin
      real(kind=wp),             allocatable :: S_svd(:,:)
      real(kind=wp),             allocatable :: Swrk(:,:)
      real(kind=wp),             allocatable :: Sigma(:)
      real(kind=wp),             allocatable :: V(:,:), W(:,:)
      class(abstract_vector),    allocatable :: Uwrk(:)

      rkc = size(X%U)
      rko = size(Y%U)
      allocate(S_svd(rkc,rko))
      ! scratch arrays
      rk = max(rkc, rko)
      rkmin = min(rkc, rko)
      allocate(Swrk(rk,rk))
      allocate(Uwrk(rk), source=T(1))

      ! compute inner product with Gramian bases
      call mat_mult(S_svd, Y%U, X%U)

      ! compute Cholesky factors of Y update LR factor with Cholesky factor
      Swrk = 0.0_wp
      call sqrtm(Swrk(1:rkc,1:rkc), Y%S)
      call mat_zero(Uwrk)
      call mat_mult(Uwrk(1:rkc), Y%U, Swrk(1:rkc,1:rkc))
      call mat_copy(Y%U, Uwrk(1:rkc))
      ! Update data matrix
      S_svd = matmul(Swrk(1:rkc,1:rkc), S_svd)

      ! compute Cholesky factors of X update LR factor with Cholesky factor
      Swrk = 0.0_wp
      call sqrtm(Swrk(1:rko,1:rko), X%S)
      call mat_zero(Uwrk)
      call mat_mult(Uwrk(1:rkc), X%U, Swrk(1:rkc,1:rkc))
      call mat_copy(X%U, Uwrk(1:rkc))
      ! Update data matrix
      S_svd = matmul(S_svd, Swrk(1:rko,1:rko))

      ! Compute BT
      allocate(V(rkc,rkc)); allocate(W(rko,rko)); allocate(Sigma(rkmin))
      call svd(S_svd, V, Sigma, W)

      ! We truncate in case come singular values are very small
      s_inv: do i = 1, rkmin
         if (Sigma(i) < atol) then
            exit s_inv
            rk = i-1
         end if 
         Sigma(i) = 1/sqrt(Sigma(i))
      enddo s_inv
      
      call mat_mult( T, X%U, matmul(          W(:, 1:rk),  diag(sigma(1:rk))))
      call mat_mult(ST, Y%U, matmul(transpose(V(1:rk, :)), diag(sigma(1:rk))))

      return
   end subroutine

   subroutine ROM_Petrov_Galerkin_Projection(romLTI, LTI, T, ST)
      !! Computes the Reduced-Order of the input LTI dynamical system via Petrov-Galerkin projection using 
      !! the biorthogonal projection bases \( \mathbf{V} \) and \( \mathbf{W} \) with 
      !! \( \mathbf{W}^T \mathbf{V} = \mathbf{I} \).
      !! 
      !! Given an LTI system defined by the matrices \( \mathbf{A}, \mathbf{B}, \mathbf{C}, \mathbf{D}\), 
      !! the matrices \( \hat{\mathbf{A}}, \hat{\mathbf{B}}, \hat{\mathbf{C}}, \hat{\mathbf{D}}\) of the 
      !! projected LTI system are given by:
      !! \[
      !!     \hat{\mathbf{A}} = \mathbf{W}^T \mathbf{A} \mathbf{V}, \qquad
      !!     \hat{\mathbf{B}} = \mathbf{W}^T \mathbf{B}, \qquad
      !!     \hat{\mathbf{C}} = \mathbf{C} \mathbf{V}, \qquad
      !!     \hat{\mathbf{D}} = \mathbf{D} .
      !! \]
      class(abstract_ROM_lti_system),   intent(out)    :: romLTI
      !! Reduced-order LTI
      class(abstract_lti_system),       intent(in)     :: LTI
      !! Large-scale LTI to project
      class(abstract_vector),           intent(inout)  :: T(:)
      !! Balancing transformation
      class(abstract_vector),           intent(in)     :: ST(:)
      !! Inverse balancing transformation

      ! internal variables
      integer                                       :: i, rk
      class(abstract_vector),           allocatable :: Uwrk(:)

      rk = size(T)
      allocate(Uwrk(rk), source=T(1)); call mat_zero(Uwrk)

      do i = 1, rk
         call LTI%A%matvec(Uwrk(i), T(i))
      end do
      call mat_mult(romLTI%A, ST, Uwrk)
      call mat_mult(romLTI%B, ST, LTI%B)
      call mat_mult(romLTI%C, LTI%CT, T)
      romLTI%D = LTI%D

   end subroutine ROM_Petrov_Galerkin_Projection

   subroutine ROM_Galerkin_Projection(romLTI, LTI, T)
      !! Computes the Reduced-Order of the input LTI dynamical system via Galerkin projection using 
      !! the orthogonal projection basis \( \mathbf{V} \) with \( \mathbf{V}^T \mathbf{V} = \mathbf{I} \).
      !! 
      !! Given an LTI system defined by the matrices \( \mathbf{A}, \mathbf{B}, \mathbf{C}, \mathbf{D}\), 
      !! the matrices \( \hat{\mathbf{A}}, \hat{\mathbf{B}}, \hat{\mathbf{C}}, \hat{\mathbf{D}}\) of the projected LTI system is given by:
      !! \[
      !!     \hat{\mathbf{A}} = \mathbf{V}^T \mathbf{A} \mathbf{V}, \qquad
      !!     \hat{\mathbf{B}} = \mathbf{V}^T \mathbf{B}, \qquad
      !!     \hat{\mathbf{C}} = \mathbf{C} \mathbf{V}, \qquad
      !!     \hat{\mathbf{D}} = \mathbf{D} .
      !! \]
      class(abstract_ROM_lti_system),   intent(out)    :: romLTI
      !! Reduced-order LTI
      class(abstract_lti_system),       intent(in)     :: LTI
      !! Large-scale LTI to project
      class(abstract_vector),           intent(inout)  :: T(:)
      !! Balancing transformation

      call ROM_Petrov_Galerkin_Projection(romLTI, LTI, T, T)

      return
   end subroutine ROM_Galerkin_Projection

end module LightROM_utils
