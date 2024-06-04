module LightROM_LyapunovUtils
   use LightKrylov
   use LightKrylov, only: wp => dp
   use LightKrylov_AbstractVectors
   implicit none

   private
   ! Matrix operations for abstract vector types

   public :: apply_outerprod

   interface apply_outerprod
      module procedure apply_outerprod_matrix_rdp
   end interface

   contains

      subroutine apply_outerprod_matrix_rdp(C,A,B)
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
         class(abstract_vector_rdp), intent(out) :: C(:)
         class(abstract_vector_rdp), intent(in)  :: A(:)
         class(abstract_vector_rdp), intent(in)  :: B(:)
         ! Intermediate basis
         real(wp),      allocatable :: wrk(:,:)

         allocate(wrk(1:size(A),1:size(B)))
         call innerprod_matrix(wrk, A, B)
         block
            class(abstract_vector_rdp), allocatable :: Xwrk(:)
            call linear_combination(Xwrk, A, wrk)
            call copy_basis(C, Xwrk)
         end block
         return
      end subroutine apply_outerprod_matrix_rdp

end module LightROM_LyapunovUtils