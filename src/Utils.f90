module LightROM_Utils
   use LightKrylov
   use LightKrylov, only : dp, wp => dp
   use LightKrylov_AbstractVectors
   use LightKrylov_Utils
   use LightROM_AbstractLTIsystems

   use stdlib_linalg, only : eye, diag, is_symmetric
   use stdlib_optval, only : optval
   implicit none 

   private
   public :: Balancing_Transformation
   public :: ROM_Petrov_Galerkin_Projection
   public :: ROM_Galerkin_Projection

   public :: zero_basis, sqrtm

   interface zero_basis
      module procedure zero_basis_rdp
   end interface

   interface sqrtm
      module procedure sqrtm_rdp
   end interface

   interface Balancing_Transformation
      module procedure Balancing_Transformation_rdp
   end interface

   interface ROM_Petrov_Galerkin_Projection
      module procedure ROM_Petrov_Galerkin_Projection_rdp
   end interface

   interface ROM_Galerkin_Projection
      module procedure ROM_Galerkin_Projection_rdp
   end interface

contains

   subroutine zero_basis_rdp(X)
      class(abstract_vector_rdp), intent(inout) :: X(:)
      integer :: i
      do i = 1, size(X)
         call X(i)%zero()
      end do

      return
   end subroutine

   subroutine sqrtm_rdp(sqrtmX, X)
      real(dp), intent(in)  :: X(:,:)
      real(dp), intent(out) :: sqrtmX(size(X,1),size(X,1))
      ! internals
      real(dp) :: lambda(size(X,1))
      real(dp) :: V(size(X,1), size(X,1))
      logical :: symmetric
      integer :: i

      ! Check if the matrix is symmetric
      if (.not. is_symmetric(X)) then
        write(*,*) "Error: Input matrix is not symmetric"
        STOP
      end if

      ! Perform eigenvalue decomposition
      call eigh(X, V, lambda)

      ! Check if the matrix is positive definite (up to tol)
      do i = 1, size(lambda)
         if (abs(lambda(i)) .gt. 1e-12_wp) then
            if (lambda(i) .gt. 0.0_wp) then
               lambda(i) = sqrt(lambda(i))
            else
               write(*,*) "Error: Input matrix is not positive definite to tolerance"
               STOP
            end if
         else
            lambda(i) = sqrt(abs(lambda(i)))
            write(*,*) "Warning: Input matrix is singular to tolerance"
         end if
      end do

      ! Reconstruct the square root matrix
      sqrtmX = matmul(V, matmul(diag(lambda), transpose(V)))  

      return
   end subroutine

   subroutine Balancing_Transformation_rdp(T, S, Tinv, Xc, Yo)
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
      class(abstract_vector_rdp),          intent(out)   :: T(:)
      !! Balancing transformation
      real(wp),                            intent(out)   :: S(:)
      !! Singular values of the BT
      class(abstract_vector_rdp),          intent(out)   :: Tinv(:)
      !! Inverse balancing transformation
      class(abstract_vector_rdp),          intent(in)    :: Xc(:)
      !! Low-rank representation of the Controllability Gramian
      class(abstract_vector_rdp),          intent(in)    :: Yo(:)
      !! Low-rank representation of the Observability Gramian

      ! internal variables
      integer                                :: i, rkc, rko, rk, rkmin
      real(wp),                  allocatable :: LRCrossGramian(:,:)
      real(wp),                  allocatable :: Swrk(:,:)
      real(wp),                  allocatable :: Sigma(:)
      real(wp),                  allocatable :: V(:,:), W(:,:)

      rkc   = size(Xc)
      rko   = size(Yo)
      rk    = max(rkc, rko)
      rkmin = min(rkc, rko) 

      ! compute inner product with Gramian bases and compte SVD
      allocate(LRCrossGramian(rkc,rko)); allocate(V(rko,rko)); allocate(W(rkc,rkc))
      call innerprod_matrix(LRCrossGramian, Xc, Yo)
      call svd(LRCrossGramian, V, S, W)

      allocate(Sigma(rkmin))
      do i = 1, rkmin
         Sigma(i) = 1/sqrt(S(i))
      enddo
      block
         class(abstract_vector_rdp), allocatable :: Xwrk(:)
         call linear_combination(Xwrk, Yo(1:rkmin), matmul(W(1:rkmin,1:rkmin), diag(Sigma)))
         call copy_basis(T(1:rkmin), Xwrk)
         call linear_combination(Xwrk, Xc(1:rkmin), matmul(V(1:rkmin,1:rkmin), diag(Sigma)))
         call copy_basis(Tinv(1:rkmin), Xwrk)
      end block

      return
   end subroutine Balancing_Transformation_rdp

   subroutine ROM_Petrov_Galerkin_Projection_rdp(Ahat, Bhat, Chat, D, LTI, T, Tinv)
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
      real(wp),            allocatable, intent(out)    :: Ahat(:, :)
      !! Reduced-order dynamics matrix.
      real(wp),            allocatable, intent(out)    :: Bhat(:, :)
      !! Reduced-order input-to-state matrix.
      real(wp),            allocatable, intent(out)    :: Chat(:, :)
      !! Reduced-order state-to-output matrix.
      real(wp),            allocatable, intent(out)    :: D(:, :)
      !! Feed-through matrix
      class(abstract_lti_system_rdp),   intent(in)     :: LTI
      !! Large-scale LTI to project
      class(abstract_vector_rdp),       intent(in)     :: T(:)
      !! Balancing transformation
      class(abstract_vector_rdp),       intent(in)     :: Tinv(:)
      !! Inverse balancing transformation

      ! internal variables
      integer                                          :: i, rk, rkc, rkb
      class(abstract_vector_rdp),       allocatable    :: Uwrk(:)
      real(wp),                         allocatable    :: Cwrk(:, :)

      rk  = size(T)
      rkb = size(LTI%B)
      rkc = size(LTI%CT)
      allocate(Uwrk(rk), source=T(1)); call zero_basis(Uwrk)
      allocate(Ahat(1:rk, 1:rk ));                  Ahat = 0.0_dp
      allocate(Bhat(1:rk, 1:rkb));                  Bhat = 0.0_wp
      allocate(Cwrk(1:rk, 1:rkc));                  Cwrk = 0.0_wp
      allocate(Chat(1:rkc,1:rk ));                  Chat = 0.0_wp
      allocate(D(1:size(LTI%D,1),1:size(LTI%D,2))); D    = 0.0_wp

      do i = 1, rk
         call LTI%A%matvec(Tinv(i), Uwrk(i))
      end do
      call innerprod_matrix(Ahat, T, Uwrk)
      call innerprod_matrix(Bhat, T, LTI%B)
      call innerprod_matrix(Cwrk, LTI%CT, Tinv)
      Chat = transpose(Cwrk)
      D = LTI%D

   end subroutine ROM_Petrov_Galerkin_Projection_rdp

   subroutine ROM_Galerkin_Projection_rdp(Ahat, Bhat, Chat, D, LTI, T)
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
      real(wp),            allocatable, intent(out)    :: Ahat(:, :)
      !! Reduced-order dynamics matrix.
      real(wp),            allocatable, intent(out)    :: Bhat(:, :)
      !! Reduced-order input-to-state matrix.
      real(wp),            allocatable, intent(out)    :: Chat(:, :)
      !! Reduced-order state-to-output matrix.
      real(wp),            allocatable, intent(out)    :: D(:, :)
      !! Feed-through matrix
      class(abstract_lti_system_rdp),   intent(in)     :: LTI
      !! Large-scale LTI to project
      class(abstract_vector_rdp),       intent(inout)  :: T(:)
      !! Balancing transformation

      call ROM_Petrov_Galerkin_Projection(Ahat, Bhat, Chat, D, LTI, T, T)

      return
   end subroutine ROM_Galerkin_Projection_rdp

end module LightROM_Utils
