module Laplacian2D_LTI_Riccati
   use Laplacian2D_LTI_Base
   !> LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov_utils
   !> Standard Libraries.
   use stdlib_linalg, only: eye
   use stdlib_optval, only: optval
   implicit none

   private
   public :: CARE

contains

   function CARE(X,A,Q,BRinvBT) result(Y)
      real(kind=wp), dimension(n,n) :: X, A, Q, BRinvBT, Y
      Y = matmul(transpose(A), X) + matmul(X, A) + Q - matmul(X, matmul(BRinvBT, X))
   end function CARE

end module Laplacian2D_LTI_Riccati