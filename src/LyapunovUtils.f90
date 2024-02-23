module LightROM_LyapunovUtils
   Use LightKrylov
   Use Lightrom_wplib
   implicit none

   private
   !> Matrix operations for abstract vector types
   public :: apply_outerproduct

   !------------------------------
   !-----     INTERFACES     -----
   !------------------------------

   contains

      subroutine apply_outerproduct(C,A,B)
         ! Compute the matrix product C = Q @ B where
         !     Q = A @ A.T is the outer product of A with itself
         ! with
         !     C: abstract vector type Krylov basis :: size nxr
         !     A: abstract vector type Krylov basis :: size nxm
         !     B: abstract vector type Krylov basis :: size nxr
         ! In order to avoid building Q (nxn), we compute sequentially
         !     C = A @ ( A.T @ B )
         class(abstract_vector) , intent(out) :: C(:)
         class(abstract_vector) , intent(in)  :: A(:)
         class(abstract_vector) , intent(in)  :: B(:)
         !> Intermediate matrix
         real(kind=wp), allocatable :: wrk(:,:)
         allocate(wrk(1:size(A),1:size(B)))
         call mat_mult(wrk, A, B)
         call mat_mult(C, A, wrk)
         deallocate(wrk)
         return
      end subroutine apply_outerproduct

end module LightROM_LyapunovUtils