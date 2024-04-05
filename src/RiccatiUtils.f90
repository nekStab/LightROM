module LightROM_RiccatiUtils
   use LightKrylov
   use LightKrylov_utils, only : assert_shape
   implicit none

   private
   public :: apply_outerproduct_w, apply_p_outerproduct_w

   !------------------------------
   !-----     INTERFACES     -----
   !------------------------------

contains

   subroutine apply_outerproduct_w(basis_out, basis_in, B, W)
      class(abstract_vector),     intent(out)  :: basis_out(:)
      class(abstract_vector),     intent(in)   :: basis_in(:)
      class(abstract_vector),     intent(in)   :: B(:)
      real(kind=wp),              intent(in)   :: W(:,:)
      ! internals
      integer                                  :: p, rk
      real(kind=wp),               allocatable :: wrk(:,:)

      p  = size(B)
      rk = size(basis_in)
      allocate(wrk(1:p,1:rk)); wrk = 0.0_wp

      call assert_shape(W, (/ p, p /), 'apply_outerproduct_w', 'W')

      call mat_zero(basis_out)
      call mat_mult(wrk, B, basis_in)
      call mat_mult(basis_out, B, matmul(W, wrk))
   
      return   
   end subroutine apply_outerproduct_w

   subroutine apply_p_outerproduct_w(mat_out, UL, UR, B, W)
      real(kind=wp),              intent(out)  :: mat_out(:,:)
      class(abstract_vector),     intent(in)   :: UL(:)
      class(abstract_vector),     intent(in)   :: UR(:)
      class(abstract_vector),     intent(in)   :: B(:)
      real(kind=wp),              intent(in)   :: W(:,:)
      ! internals
      real(kind=wp)                            :: BTUR(size(B),size(UR))
      real(kind=wp)                            :: ULTB(size(UL),size(B))

      call assert_shape(mat_out, (/ size(UL), size(UR) /), 'apply_p_outerproduct_w', 'mat_out')
      
      BTUR = 0.0_wp; ULTB = 0.0_wp; mat_out = 0.0_wp
      call mat_mult(BTUR, B, UR)
      call mat_mult(ULTB, UL, B)
      
      mat_out = matmul( ULTB, matmul( W, BTUR ) )
   
      return   
   end subroutine apply_p_outerproduct_w

end module LightROM_RiccatiUtils