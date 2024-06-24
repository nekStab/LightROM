module LightROM_LyapunovUtils
   ! LightKrylov
   use LightKrylov
   use LightKrylov, only: wp => dp
   use LightKrylov_Constants
   use LightKrylov_AbstractVectors
   ! LightROM
   use LightROM_AbstractLTIsystems
   
   implicit none

   ! module name
   private :: this_module
   character*128, parameter :: this_module = 'LightROM_LyapunovUtils'

   ! Matrix operations for abstract vector types
   public :: apply_outerprod

   interface apply_outerprod
      module procedure apply_outerprod_vector_rdp
      module procedure apply_outerprod_basis_rdp
   end interface

contains

   subroutine apply_outerprod_vector_rdp(c,A,b)
      !! Computes the matrix product \( \mathbf{c} = \mathbf{Q} \mathbf{b} \) where
      !! \( \mathbf{Q} = \mathbf{A} \mathbf{A}^T \) is the outer product of \( \mathbf{A} \)
      !! with itself with
      !!
      !! - \( \mathbf{C} \): `abstract vector`
      !! - \( \mathbf{A} \): `abstract vector` type Krylov basis of size (n x m)
      !! - \( \mathbf{B} \): `abstract vector`
      !!
      !! In order to avoid building \( \mathbf{Q} \) (n x n), we compute sequentially
      !! \( \mathbf{c} = \mathbf{A} ( \mathbf{A}^T \mathbf{b} ) \)
      class(abstract_vector_rdp), intent(out) :: c
      class(abstract_vector_rdp), intent(in)  :: A(:)
      class(abstract_vector_rdp), intent(in)  :: b
      ! Intermediate basis
      real(wp) :: wrk(size(A))
      wrk = zero_rdp
      call innerprod(wrk, A, b)
      block
         class(abstract_vector_rdp), allocatable :: xwrk
         call linear_combination(xwrk, A, wrk)
         call c%zero(); call c%add(xwrk)
      end block
      return
   end subroutine apply_outerprod_vector_rdp

   subroutine apply_outerprod_basis_rdp(C,A,B)
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
      real(wp) :: wrk(size(A),size(B))
      wrk = zero_rdp
      call innerprod(wrk, A, B)
      block
         class(abstract_vector_rdp), allocatable :: Xwrk(:)
         call linear_combination(Xwrk, A, wrk)
         call copy_basis(C, Xwrk)
      end block
      return
   end subroutine apply_outerprod_basis_rdp

end module LightROM_LyapunovUtils