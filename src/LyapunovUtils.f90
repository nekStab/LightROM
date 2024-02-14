module LightROM_LyapunovUtils
   Use LightKrylov
   implicit none

   private
   !> Matrix products for abstract vector types
   public :: mat_mult, apply_rhs

   !------------------------------
   !-----     INTERFACES     -----
   !------------------------------

   interface mat_mult
      module procedure mat_mult_direct
      module procedure mat_mult_transpose
   end interface mat_mult

   contains

      !-----------------------------------
      !-----                         -----
      !-----     MATRIX PRODUCTS     -----
      !-----                         -----
      !-----------------------------------
   
      subroutine mat_mult_direct(C,A,B)
        ! Compute the matrix product C = A @ B with
        !     C: abstract vector type Krylov basis :: size nxq
        !     A: abstract vector type Krylov basis :: size nxr
        !     B: real matrix                       :: size rxq
        class(abstract_vector) , intent(out) :: C(:)   ! result
        class(abstract_vector) , intent(in)  :: A(:)   ! krylov basis
        real(kind=wp)          , intent(in)  :: B(:,:) ! coefficient matrix
        !> Local variables
        class(abstract_vector) , allocatable :: wrk
        integer :: i
      
        !> Check sizes.
        if (size(A) .ne. size(B,1)) then
           write(*,*) "INFO : Abstract vector basis A and coefficient matrix B have incompatible sizes for the product A @ B."
           return
        endif
        allocate(wrk, source=A(1))
        !> Compute product column-wise
        do i = 1, size(B,2)
           call get_vec(wrk, A, B(:, i))
           call C(i)%axpby(0.0_wp, wrk, 1.0_wp)
        enddo
        deallocate(wrk)
        return
      end subroutine mat_mult_direct

      subroutine mat_mult_transpose(C,A,B)
         ! Compute the matrix product C = A.T @ B with
         !     C: real matrix                       :: size rxq
         !     A: abstract vector type Krylov basis :: size nxr
         !     B: abstract vector type Krylov basis :: size nxq
         class(abstract_vector) , intent(in)  :: A(:)   ! krylov basis
         class(abstract_vector) , intent(in)  :: B(:)   ! krylov basis
         real(kind=wp)          , intent(out) :: C(size(A),size(B)) ! result
         !> Local variables
         integer :: i, j

         !> Compute product column-wise
         C = 0.0_wp
         do i = 1, size(A)
            do j = 1, size(B)
               C(i,j) = A(i)%dot(B(j))
            enddo
         enddo
         return

      end subroutine mat_mult_transpose

      subroutine apply_rhs(C,A,B)
         ! Compute the matrix product C = Q @ B where
         !     Q = A @ A.T 
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
         call mat_mult_transpose(wrk, A, B)
         call mat_mult_direct(C, A, wrk)
         deallocate(wrk)
         return
      end subroutine apply_rhs

end module LightROM_LyapunovUtils