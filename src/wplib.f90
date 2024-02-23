module LightROM_wplib
   Use LightKrylov
   Use lightkrylov_utils

   !> Fortran standard library.
   use stdlib_optval, only: optval

   implicit none

   private
   !> Matrix operations for abstract vector types
   public :: kexpm, expm, factorial

contains

   subroutine kexpm(C, A, B, tau, tol, info, nkryl)
      !> Krylov basis of the approximate value of expm(tA) @ B
      class(abstract_vector), intent(out) :: C(:)
      !> Linear operator to be exponentiated
      class(abstract_linop),  intent(in) :: A
      !> Krylov basis on which to apply expm(tA)
      class(abstract_vector), intent(in) :: B(:)
      !> time horizon for exponentiation
      real(kind=wp),          intent(in) :: tau
      !> solution tolerance based on error estimates
      real(kind=wp),          intent(in) :: tol
      !> Information flag
      integer,                intent(out) :: info
      !> Optional
      integer, optional,      intent(in)  :: nkryl
      integer                             :: nk
      
      !> internals
      integer, parameter :: kmax = 10
      integer :: k, p, kpm, kp, kpp
      !> Arnoldi factorisation
      class(abstract_vector), allocatable :: X(:)
      real(kind=wp), allocatable          :: H(:,:)
      !> Normalisation & temp arrays
      real(kind=wp), allocatable          :: R(:,:)
      real(kind=wp), allocatable          :: E(:,:)
      class(abstract_vector), allocatable :: wrk(:)

      !> determine block size
      p = size(B)

      !> Optional arguemnts
      nk = optval(nkryl, kmax)

      ! allocate memory
      allocate(R(1:p,1:p)); allocate(H(1:p*(nk+1),1:p*nk));
      allocate(E(1:nk,1:nk));
      allocate(wrk(1:nk), source=B(1))
      allocate(X(1:nk+1), source=B(1))

      !> normalize input vector and set initialise Krylov subspace
      R = 0.0_wp
      call mat_zero(wrk)
      call mat_copy(wrk(1:p),B(1:p))
      call qr_factorization(wrk(1:p),R(1:p,1:p),info)
      call initialize_krylov_subspace(X,wrk(1:p))
      H = 0.0_wp

      do k = 1, nk
         kpm = (k-1)*p
         kp  = kpm + p
         kpp = kp  + p
         !> reset wrk arrays
         E = 0.0_wp; call mat_zero(wrk)
         !> compute kth stop of the Arnoldi factorization
         call arnoldi_factorization(A, X(1:kpp), H(1:kpp,1:kp), info, kstart=k, kend=k, block_size=p)
         !> compute the (dense) matrix exponential of the Hessenberg matrix
         call expm(E(1:kp,1:kp),tau*H(1:kp,1:kp))
         !> project back into original space
         call mat_mult(wrk(1:kp),X(1:kp),E(1:kp,1:kp))
         !> extract relevant part and scale
         call mat_mult(C(1:p),wrk(1:p),R(1:p,1:p))
      end do
      
   end subroutine kexpm
   
   subroutine expm(E,A)
   
   !*****************************************************************************80
   !
   !! EXPM_r is essentially MATLAB's built-in matrix exponential algorithm.
   !
   !  Licensing:
   !
   !    This code is distributed under the MIT license.
   !
   !  Modified:
   !
   !    26 November 2011
   !
   !  Author:
   !
   !    Cleve Moler, Charles Van Loan
   !
   !  Reference:
   !
   !    Cleve Moler, Charles VanLoan,
   !    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
   !    Twenty-Five Years Later,
   !    SIAM Review,
   !    Volume 45, Number 1, March 2003, pages 3-49.
   !
   !  Input:
   !
   !    integer N, the dimension of the matrix.
   !
   !    real(kind=rk) A(N,N), the matrix.
   !
   !  Output:
   !
   !    real(kind=rk) E(N,N), the estimate for exp(A).
   !
     implicit none
   
     !> arguments
     real(kind=wp), intent(in)  :: a(:,:)
     real(kind=wp), intent(out) :: e(:,:)
     !> Internal
     real(kind=wp), allocatable :: A2(:,:), D(:,:), X(:,:), invD(:,:)
     real(kind=wp)              :: a_norm, c, norm_li
     integer, parameter         :: q = 6
     integer                    :: n, ee, k, s
     logical                    :: p

     n = size(A,1)
     allocate(A2(1:n,1:n)); allocate(X(1:n,1:n))
     allocate(D(1:n,1:n));  allocate(invD(1:n,1:n))

     !> Make a copy of the matrix.
     a2(1:n,1:n) = a(1:n,1:n)
     !> Compute the L-infinity norm.
     a_norm = norm_linf(A2)
     !> Determine a scaling factor for the matrix.
     ee = int(log2(a_norm)) + 1
     s  = max(0, ee + 1)

     a2 = a2 / 2.0_wp**s
     x  = a2
   
     c = 0.5_wp
   
     e = 0.0_wp; forall (k=1:n) e(k, k) = 1.0_wp
     e = e + c * a2
   
     d = 0.0_wp; forall (k=1:n) d(k, k) = 1.0_wp
     d = d - c * a2
   
     p = .true.
   
     do k = 2, q
       c = c*( q - k + 1)/(k * ( 2 * q - k + 1 ))
       x = matmul( a2, x )
       e = e + c * x
       if ( p ) then
         d = d + c * x
       else
         d = d - c * x
       end if
       p = .not. p
     end do
   !
   !  E -> inverse(D) * E
   !
     invD = D
     call inv(invD)
     E = matmul(invD,E)
   !
   !  E -> E^(2*S)
   !
     do k = 1, s
       e = matmul ( e, e )
     end do

     deallocate(A2); deallocate(X); deallocate(D); deallocate(invD)
   
     return
   end subroutine expm

   function log2(x) result(y)
     implicit none
     real(kind=wp), intent(in) :: x
     real(kind=wp) :: y
     y = log(x) / log(2.0_wp)
   end function log2

   function norm_linf(A) result(norm)
      implicit none   
      real(kind=wp), intent(in) :: A(:,:)
      !> internals
      integer       :: i, n
      real(kind=wp) :: row_sum, norm
      !> initialize
      norm = 0.0_wp
      n = size(A,1)
      do i = 1, n
      row_sum = sum ( abs ( A(i,1:n) ) )
      norm = max ( norm, row_sum )
      end do
   end function norm_linf

   function factorial(n) result(y)
      implicit none
      integer, intent(in) :: n 
      integer :: i, y
      y = 1
      do i = 2, n
         y = y*i 
      end do
   end function factorial

end module LightROM_wplib