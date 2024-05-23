module LightROM_utils
   use LightKrylov
   use LightKrylov_utils
   use LightROM_AbstractLTIsystems

   use stdlib_linalg, only : eye, diag
   use stdlib_optval, only : optval
   implicit none 

   private
   public Balancing_Transformation, ROM_Petrov_Galerkin_Projection, ROM_Galerkin_Projection

contains

   subroutine Balancing_Transformation(T,S,Tinv,Xc,Yo)
      !! Computes the the biorthogonal balancing transformation \( \mathbf{T}, \mathbf{T}^{-1} \) from the
      !! low-rank representation of the SVD of the controllability and observability Gramians, \( \mathbf{W}_c \) 
      !! and \( \mathbf{W}_o \) respectively, given as:
      !! \[ \begin{align}
      !!    \mathbf{W}_c &= \mathbf{X}_c \mathbf{X}_c^T \\
      !!    \mathbf{W}_o &= \mathbf{Y}_o \mathbf{Y}_o^T
      !! \end{align} \]
      !!
      !! Given the SVD of the cross-Gramian:
      !! $$ \mathbf{X}_c^T \mathbf{Y}_o = \mathbf{U} \mathbf{S} \mathbf{V}^T $$
      !! the balancing transformation and its inverse are given by:
      !! \[ \begin{align}
      !!            \mathbf{T}   &= \mathbf{X}_o \mathbf{S}_o^{1/2} \mathbf{V} \mathbf{S}^{-1/2} \\
      !!            \mathbf{Tinv}^T &= \mathbf{Y}_c \mathbf{S}_c^{1/2} \mathbf{U} \mathbf{S}^{-1/2} 
      !! \end{align} \]
      !! Note: In the current implementation, the numerical rank of the SVD is not considered.
      class(abstract_vector),              intent(out)   :: T(:)
      !! Balancing transformation
      real(kind=wp),                       intent(out)   :: S(:)
      !! Singular values of the BT
      class(abstract_vector),              intent(out)   :: Tinv(:)
      !! Inverse balancing transformation
      class(abstract_vector),              intent(in)    :: Xc(:)
      !! Low-rank representation of the Controllability Gramian
      class(abstract_vector),              intent(in)    :: Yo(:)
      !! Low-rank representation of the Observability Gramian

      ! internal variables
      integer                                :: i, rkc, rko, rk, rkmin
      real(kind=wp),             allocatable :: LRCrossGramian(:,:)
      real(kind=wp),             allocatable :: Swrk(:,:)
      real(kind=wp),             allocatable :: Sigma(:)
      real(kind=wp),             allocatable :: V(:,:), W(:,:)

      rkc   = size(Xc)
      rko   = size(Yo)
      rk    = max(rkc, rko)
      rkmin = min(rkc, rko) 

      ! compute inner product with Gramian bases and compte SVD
      allocate(LRCrossGramian(rkc,rko)); allocate(V(rko,rko)); allocate(W(rkc,rkc))
      call mat_mult(LRCrossGramian, Xc, Yo)
      call svd(LRCrossGramian, V, S, W)

      allocate(Sigma(rkmin))
      do i = 1, rkmin
         Sigma(i) = 1/sqrt(S(i))
      enddo
         
      call mat_mult(T(1:rkmin),    Yo(1:rkmin), matmul(W(1:rkmin,1:rkmin), diag(Sigma)))
      call mat_mult(Tinv(1:rkmin), Xc(1:rkmin), matmul(V(1:rkmin,1:rkmin), diag(Sigma)))

      return
   end subroutine Balancing_Transformation

   subroutine ROM_Petrov_Galerkin_Projection(Ahat, Bhat, Chat, D, LTI, T, Tinv)
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
      real(kind=wp),       allocatable, intent(out)    :: Ahat(:, :)
      !! Reduced-order dynamics matrix.
      real(kind=wp),       allocatable, intent(out)    :: Bhat(:, :)
      !! Reduced-order input-to-state matrix.
      real(kind=wp),       allocatable, intent(out)    :: Chat(:, :)
      !! Reduced-order state-to-output matrix.
      real(kind=wp),       allocatable, intent(out)    :: D(:, :)
      !! Feed-through matrix
      class(abstract_lti_system),       intent(in)     :: LTI
      !! Large-scale LTI to project
      class(abstract_vector),           intent(in)     :: T(:)
      !! Balancing transformation
      class(abstract_vector),           intent(in)     :: Tinv(:)
      !! Inverse balancing transformation

      ! internal variables
      integer                                          :: i, rk, rkc, rkb
      class(abstract_vector),           allocatable    :: Uwrk(:)
      real(kind=wp),                    allocatable    :: Cwrk(:, :)

      rk  = size(T)
      rkb = size(LTI%B)
      rkc = size(LTI%CT)
      allocate(Uwrk(rk), source=T(1)); call mat_zero(Uwrk)
      allocate(Ahat(1:rk, 1:rk ));                  Ahat = 0.0_wp
      allocate(Bhat(1:rk, 1:rkb));                  Bhat = 0.0_wp
      allocate(Cwrk(1:rk, 1:rkc));                  Cwrk = 0.0_wp
      allocate(Chat(1:rkc,1:rk ));                  Chat = 0.0_wp
      allocate(D(1:size(LTI%D,1),1:size(LTI%D,2))); D    = 0.0_wp

      do i = 1, rk
         call LTI%A%matvec(Tinv(i), Uwrk(i))
      end do
      call mat_mult(Ahat, T, Uwrk)
      call mat_mult(Bhat, T, LTI%B)
      call mat_mult(Cwrk, LTI%CT, Tinv)
      Chat = transpose(Cwrk)
      D = LTI%D

   end subroutine ROM_Petrov_Galerkin_Projection

   subroutine ROM_Galerkin_Projection(Ahat, Bhat, Chat, D, LTI, T)
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
      real(kind=wp),       allocatable, intent(out)    :: Ahat(:, :)
      !! Reduced-order dynamics matrix.
      real(kind=wp),       allocatable, intent(out)    :: Bhat(:, :)
      !! Reduced-order input-to-state matrix.
      real(kind=wp),       allocatable, intent(out)    :: Chat(:, :)
      !! Reduced-order state-to-output matrix.
      real(kind=wp),       allocatable, intent(out)    :: D(:, :)
      !! Feed-through matrix
      class(abstract_lti_system),       intent(in)     :: LTI
      !! Large-scale LTI to project
      class(abstract_vector),           intent(inout)  :: T(:)
      !! Balancing transformation

      call ROM_Petrov_Galerkin_Projection(Ahat, Bhat, Chat, D, LTI, T, T)

      return
   end subroutine ROM_Galerkin_Projection

end module LightROM_utils
