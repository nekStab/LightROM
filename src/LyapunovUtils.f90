module LightROM_LyapunovUtils
   Use LightKrylov
   Use LightKrylov_expmlib
   implicit none

   private
   ! Matrix operations for abstract vector types
   public :: apply_outerproduct

   contains

      subroutine apply_outerproduct(C,A,B)
         !! Computes the matrix product \( \mathbf{C} = \mathbf{Q} \mathbf{B} \) where
         !! \( \mathbf{Q} = \mathbf{A} \mathbf{A}^T \) is the outer product of \( \mathbf{A} \)
         !! with itself with
         !!
         !! - \( \mathbf{C} \): `abstract vector` type Krylov basis of size (n x r)
         !! - \( \mathbf{A} \): `abstract vector` type Krylov basis of size (n x m)
         !! - \( \mathbf{B} \): `abstract vector` type Krylov basis of size (n x r)
         !!
         !! In order to avoid building \( \mathbf{Q} \) (n x n), we compute sequentially
         !! \( \mathbf{C} = \mathbf{A} ( \mathbf{A}^T \mathbf{B} ) \)
         class(abstract_vector) , intent(out) :: C(:)
         class(abstract_vector) , intent(in)  :: A(:)
         class(abstract_vector) , intent(in)  :: B(:)
         ! Intermediate basis
         real(kind=wp), allocatable :: wrk(:,:)
         allocate(wrk(1:size(A),1:size(B)))
         call mat_mult(wrk, A, B)
         call mat_mult(C, A, wrk)
         deallocate(wrk)
         return
      end subroutine apply_outerproduct

end module LightROM_LyapunovUtils