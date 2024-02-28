module LightROM_expmlib
   Use LightKrylov
   Use lightkrylov_utils

   !> Fortran standard library.
   use stdlib_optval, only: optval

   implicit none

   private
   !> Matrix operations for abstract vector types
   public :: kexpm, expm, norm_linf, norm_fro, factorial

   interface kexpm
      module procedure kexpm_vec
      module procedure kexpm_mat
   end interface

contains

   subroutine kexpm_vec(c, A, b, tau, tol, info, verbosity, nkryl)

   !=======================================================================================
   ! Krylov-based approximation of the action the exponential propagator on a matrix B
   !=======================================================================================
   !
   !  Approximate v = expm(tau * A) @ b using Krylov subspace techniques
   !
   !=======================================================================================
      !> Krylov basis of the approximate value of expm(tA) @ B
      class(abstract_vector), intent(out) :: c
      !> Linear operator to be exponentiated
      class(abstract_linop),  intent(in) :: A
      !> Krylov basis on which to apply expm(tA)
      class(abstract_vector), intent(in) :: b
      !> time horizon for exponentiation
      real(kind=wp),          intent(in) :: tau
      !> solution tolerance based on error estimates
      real(kind=wp),          intent(in) :: tol
      !> Information flag
      integer,                intent(out) :: info
      !> Optional
      logical, optional,      intent(in)  :: verbosity
      logical                             :: verbose
      integer, optional,      intent(in)  :: nkryl
      integer                             :: nstep
      
      !> internals
      integer, parameter :: kmax = 20
      integer :: i, k, p, km, kp, nk
      !> Arnoldi factorisation
      class(abstract_vector), allocatable :: X(:)
      real(kind=wp), allocatable          :: H(:,:)
      !> Normalisation & temp arrays
      real(kind=wp), allocatable          :: E(:,:), invH(:,:), phi(:,:), wrk(:,:)
      real(kind=wp), allocatable          :: Id(:,:)
      class(abstract_vector), allocatable :: xwrk
      real(kind=wp)                       :: err_est, beta

      !> Optional arguemnts
      verbose = optval(verbosity, .false.)
      nstep   = optval(nkryl, kmax)
      nk      = nstep

      info = 0

      ! allocate memory
      allocate(X(1:nk+1), source=b); allocate(H(1:nk+1,1:nk))      ! Arnoldi factorization
      allocate(E(1:nk+1,1:nk))                                     ! Dense matrix exponential
      allocate(invH(1:nk,1:nk)); allocate(phi(1:nk,1:nk));         ! Correction
      allocate(Id(1:nk,1:nk)); Id = 0.0_wp; forall (i=1:nk) Id(i, i) = 1.0_wp
      ! scratch arrays
      allocate(wrk(1:nk,1:nk))
      allocate(xwrk, source=b)

      !> normalize input vector and set initialise Krylov subspace
      beta = b%norm()
      !call initialize_krylov_subspace(X,b)
      call mat_zero(X)
      call X(1)%axpby(0.0_wp, b, 1.0_wp)
      call X(1)%scal(1.0/beta)
      H = 0.0_wp

      expm_arnoldi: do k = 1, nk 
         km = k - 1
         kp = k + 1
         !> reset wrk arrays
         E = 0.0_wp
         !> compute kth stop of the Arnoldi factorization
         call arnoldi_factorization(A, X(1:kp), H(1:kp,1:k), info, kstart=k, kend=k)
         !> compute the (dense) matrix exponential of the Hessenberg matrix
         call expm(E(1:k,1:k), tau*H(1:k,1:k))
         !> compute correction based on last row of Hessenberg matrix
         wrk = 0.0_wp; invH = 0.0_wp
         invH(1:k,1:k) = tau*H(1:k,1:k)
         call inv(invH(1:k,1:k))
         wrk = E(1:k,1:k) - Id(1:k,1:k)
         phi(1:k,1:k) = matmul(invH(1:k,1:k), wrk(1:k,1:k))    ! phi(tH)
         !> add correction row(s) to exptH
         E(kp,1:k) = phi(k,1:k)
         !> project back into original space
         call get_vec(xwrk, X(1:k), E(1:k,1))
         call c%axpby(0.0_wp, xwrk, beta)
         !> cheap error estimate
         err_est = H(kp,k) * abs(phi(k,1) * beta)   ! h_{m+1,m} @ | e_mT @ phi(tH) @ e_1 @ beta |
         !> correction
         !call c%axpby(1.0_wp, X(kp), E(kp,1))
         if (err_est .lt. tol) then
            if (verbose) then
               write(*, *) 'Block Arnoldi approxmation of the action of the exp. propagator converged'
               write(*, *) '    n° of vectors:', kp
               write(*, *) '    estimated error:   ||err_est||_2 = ', err_est
               write(*, *) '    desired tolerance:           tol = ', tol
               write(*, *)
            endif
            ! --> Exit the Arnoldi iteration.
            exit expm_arnoldi
         endif
      end do expm_arnoldi

      if (err_est .gt. tol) then
         info = -1
         if (verbose) then
            write(*, *) 'Arnoldi-based approxmation of the exp. propagator did not converge'
            write(*, *) '    maximum n° of vectors reached: ', nk + 1
            write(*, *) '    estimated error:   ||err_est||_2 = ', err_est
            write(*, *) '    desired tolerance:           tol = ', tol
            write(*, *)
         endif
      endif

      return
   end subroutine kexpm_vec

   subroutine kexpm_mat(C, A, B, tau, tol, info, verbosity, nkryl)

      !=======================================================================================
      ! Krylov-based approximation of the action the exponential propagator on a matrix B
      !=======================================================================================
      !
      ! Purpose:
      ! --------
      ! Approximates the action of the exponential propagator (matrix exponential) of a linear 
      ! operator A on a given matrix B by computing the action of the exponential propagator
      ! on the projection of the operator on a (small) Krylov subspace.
      !
      ! Mathematical Formulation:
      ! -------------------------
      ! Given a linear operator A \in R^(n x n), we compute:
      ! 
      !              expm(tau*A) @ B ~ Xm_ @ expm(tau*Hm_) @ e_1 @ R
      !
      ! with:
      ! - Xm_ := [ X(:,1:m), X(:,m+1) ] is the basis generated by the m-step arnoldi process
      !       s. t. A @ X(:, m) = X(:, 1:m) @ H(1:m, m) + h_{m+1,m} * X(:, m+1).
      ! - Hm_ := [ Hm,         , 0 ]
      !          [ c @ phi(Hm) , 1 ]
      !      where 
      !           1. Hm = H(1:m, m)
      !           2. c = h_{m+1,m}*e_m.T
      !           3. phi(Z) = Z^(-1) @ (exp(Z) - I)
      ! - B := X(:,1) @ R is the QR factorization of the input matrix
      !
      ! Algorithmic Features:
      ! ---------------------
      ! - Very good approximations of the action of the exponential propagator are obtained 
      !   with very few arnoldi steps. O(10) steps are sufficient for working precision for 
      !   "small" matrices A (or tau << 1)
      ! - The small dense problem of the exponential of the Hessenberg matrix is solved using
      !   the scaling and squaring approach combined with a rational Pade approximation.
      ! - The correction to the naive Krylov approach proposed by Gallpoulos & Saad making
      !   use of the additional information in the last row of the Hessenberg matrix implemented
      !
      ! Advantages:
      ! -----------
      ! - Very accurate for "small" matrices (in terms of their norm), i.e. for small tau.
      ! - A fairly accuracte error estimate is computed based on the Hessenberg matrix to
      !   terminate iteration when needed (see Er3, p 223, Y. Saad, 1992)
      ! - Block arnoldi method allows for the right hand side to contain more than 1 vector
      !
      ! Limitations:
      ! ------------
      ! - No general computation of expm(tau*A) but action of the propagator on a specific B.
      !   The process must be restarted for each new B.
      ! - The error estimation is quite crude and conservative. In practise, the approximation
      !   is better than the estimator suggests
      !
      ! Input/Output Parameters:
      ! ------------------------
      ! C      : Approximation of the action of expm(tau*A), class(abstract_vector),   [Output]
      ! A      : Linear Operator, class(abstract_linop)                                [Input]
      ! B      : Input matrix, class(abstract_vector),                                 [Input]
      ! tau    : Integration time for the exponential propagator, real_wp              [Input]
      ! tol    : tolerance for the error estimator (not implemented yet), real_wp      [Input]
      ! info   : Exit information flag, integer,                                       [Output]
      ! nkryl  : Optional, maximum number of arnoldi steps, integer                    [Input]
      !
      !  References:
      ! -----------
      !  - Y. Saad, "Analysis of Some Krylov Subspace Approximations to the Matrix Exponent", 
      !    SIAM Journal on Numerical Analysis, Volume 29, Number 1, Feb. 1992, pages 209-228.
      !  - E. Gallopoulos & Y. Saad, "Efficient soltuion of parabolic equations by Krylov
      !    approximation methods", SIAM Journal on Scientific and Statistical Computing, 
      !    Volume 13, Number 5, September 1992, pages 1236-1264.
      !  - C. Moler & C. VanLoan, "Nineteen Dubious Ways to Compute the Exponential of a 
      !    Matrix, Twenty-Five Years Later", SIAM Review, Volume 45, Number 1, March 2003, 
      !    pages 3-49.
      !
      !=======================================================================================
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
         logical, optional,      intent(in)  :: verbosity
         logical                             :: verbose
         integer, optional,      intent(in)  :: nkryl
         integer                             :: nstep
         
         !> internals
         integer, parameter :: kmax = 20
         integer :: i, k, p, kpm, kp, kpp, nk
         !> Arnoldi factorisation
         class(abstract_vector), allocatable :: X(:)
         real(kind=wp), allocatable          :: H(:,:)
         !> Normalisation & temp arrays
         real(kind=wp), allocatable          :: R(:,:), E(:,:), invH(:,:), phi(:,:), wrk(:,:)
         real(kind=wp), allocatable          :: Id(:,:), em(:,:), ym(:,:)
         class(abstract_vector), allocatable :: Xwrk(:), Cwrk(:)
         real(kind=wp) :: err_est
   
         !> determine block size
         p = size(B)
   
         !> Optional arguemnts
         verbose = optval(verbosity, .false.)
         nstep   = optval(nkryl, kmax)
         nk      = nstep*p
   
         info = 0
   
         ! allocate memory
         allocate(R(1:p,1:p))                                              ! QR factorization
         allocate(X(1:nk+p), source=B(1)); allocate(H(1:p*(nk+1),1:p*nk))  ! Arnoldi factorization
         allocate(E(1:p*(nk+1),1:nk))                                      ! Dense matrix exponential
         allocate(invH(1:nk,1:nk)); allocate(phi(1:nk,1:nk));              ! Correction
         allocate(Id(1:nk,1:nk)); Id = 0.0_wp; forall (i=1:nk) Id(i, i) = 1.0_wp
         allocate(em(1:p,1:p)); allocate(ym(1:p,1:p))                     ! error estimate
         ! scratch arrays
         allocate(wrk(1:nk,1:nk))
         allocate(Xwrk(1:p), source=B(1)); allocate(Cwrk(1:p), source=B(1))
   
         !> normalize input vector and set initialise Krylov subspace
         R = 0.0_wp
         call mat_zero(Xwrk)
         call mat_copy(Xwrk(1:p),B(1:p))
         call qr_factorization(Xwrk(1:p),R(1:p,1:p),info)
         call initialize_krylov_subspace(X,Xwrk(1:p))
         H = 0.0_wp
   
         expm_arnoldi: do k = 1, nk 
            kpm = (k-1)*p
            kp  = kpm + p
            kpp = kp  + p
            !> reset wrk arrays
            E = 0.0_wp; call mat_zero(Xwrk)
            !> compute kth stop of the Arnoldi factorization
            call arnoldi_factorization(A, X(1:kpp), H(1:kpp,1:kp), info, kstart=k, kend=k, block_size=p)
            !> compute the (dense) matrix exponential of the Hessenberg matrix
            call expm(E(1:kp,1:kp),tau*H(1:kp,1:kp))
            !> compute correction based on last row of Hessenberg matrix
            wrk = 0.0_wp; invH = 0.0_wp
            invH(1:kp,1:kp) = tau*H(1:kp,1:kp)
            call inv(invH(1:kp,1:kp))
            wrk = E(1:kp,1:kp) - Id(1:kp,1:kp)
            phi(1:kp,1:kp) = matmul(invH(1:kp,1:kp),wrk(1:kp,1:kp))    ! phi(tH)
            !> add correction row(s) to exptH
            E(kp+1:kpp,1:kp) = phi(kpm+1:kp,1:kp)
            !> project back into original space
            call mat_zero(C)
            call mat_mult(Xwrk(1:p),X(1:kp),E(1:kp,1:p))
            call mat_mult(C(1:p),Xwrk(1:p),R(1:p,1:p))
            !> cheap error estimate
            err_est = 0.0_wp; wrk = 0.0_wp
            wrk(1:p,1:p) = matmul(phi(kpm+1:kp,1:p),R(1:p,1:p))        ! e_mT @ phi(tH) @ e_1 @ R
            em = matmul(H(kp+1:kpp,kpm+1:kp),wrk(1:p,1:p))             ! h_{m+1,m} @ wrk
            err_est = norm_fro(em)
            !> correction
            call mat_zero(Cwrk)
            call mat_mult(Cwrk(1:p),X(kp+1:kpp),E(kp+1:kpp,1:p))
            call mat_axpby(C,1.0_wp,Cwrk,1.0_wp)
            if (err_est .lt. tol) then
               if (verbose) then
                  if (p.eq.1) then
                     write(*, *) 'Arnoldi approxmation of the action of the exp. propagator converged'
                  else
                     write(*, *) 'Block Arnoldi approxmation of the action of the exp. propagator converged'
                  endif 
                  write(*, *) '    n° of vectors:', k+1, 'per input vector, total:', kpp
                  write(*, *) '    estimated error:   ||err_est||_2 = ', err_est
                  write(*, *) '    desired tolerance:           tol = ', tol
                  write(*, *)
               endif
               ! --> Exit the Arnoldi iteration.
               exit expm_arnoldi
            endif
         end do expm_arnoldi
   
         if (err_est .gt. tol) then
            info = -1
            if (verbose) then
               write(*, *) 'Arnoldi-based approxmation of the exp. propagator did not converge'
               write(*, *) '    maximum n° of vectors reached: ', nstep+1,&
                           & 'per input vector, total:', kpp
               write(*, *) '    estimated error:   ||err_est||_2 = ', err_est
               write(*, *) '    desired tolerance:           tol = ', tol
               write(*, *)
            endif
         endif
   
         return
      end subroutine kexpm_mat
   
   subroutine expm(E,A, order)
   
   !=======================================================================================
   ! Dense matrix exponential for n x n matrices
   !=======================================================================================
   !
   ! Purpose:
   ! --------
   ! Computes the matrix exponential E = expm(A) for a real-valued dense square matrix of 
   ! order n using the scaling and squaring approach, where the scaled problem is computed
   ! using a rational matrix Pade approximation of the form
   ! 
   !         R_pq(x) = [ Q_q(x) ]^(-1) @ P_p(x) 
   !
   ! where p and q are the respective orders of the polynomials P(x) and Q(x).
   !
   ! Algorithmic Features:
   ! ---------------------
   ! - Only diagonal approximations (p = q) are used as these are more efficient to compute
   !   than off-diagonal approximations for the same approximation order
   ! - 12th order diagonal Pade approximation by default
   !
   ! Advantages:
   ! -----------
   ! - Very accurate for "small" matrices (in terms of their norm). When considering the
   !   temporal evolution of an LTI system, i.e. expm(tA), the approximation is very good
   !   for small t.
   ! - Avoids numerical problems linked to the "hump" appearing for non-normal matrices A
   !
   ! Limitations:
   ! ------------
   ! - No error estimation/check based on the approximation order and the matrix norm
   !
   ! Input/Output Parameters:
   ! ------------------------
   ! E          : Approximation of the matrix exponential   [Output]
   ! A          : Real-valued matrix to exponentiate        [Input]
   ! order      : Optional, order of the polynomials in the Pade approximant. The final 
   !              approximation has the order 2*order
   !
   !  expm is an adaptation of the function R8mat_epxm1 by C. Moder & C. vanLoan
   !
   !  Original licensing:
   !    This code is distributed under the MIT license.
   !
   !  References:
   ! -----------
   !  - C. Moler & C. VanLoan, "Nineteen Dubious Ways to Compute the Exponential of a 
   !    Matrix, Twenty-Five Years Later", SIAM Review, Volume 45, Number 1, March 2003, 
   !    pages 3-49.
   !
   !=======================================================================================
     implicit none
   
     !> arguments
     real(kind=wp), intent(in)  :: A(:,:)
     real(kind=wp), intent(out) :: E(:,:)
     !> Optional
     integer, optional,      intent(in)  :: order
     integer                             :: p_order
     !> Internal
     real(kind=wp), allocatable :: A2(:,:), Q(:,:), X(:,:), invQ(:,:)
     real(kind=wp)              :: a_norm, c, norm_li
     integer                    :: n, ee, k, s
     logical                    :: p

     !> Optional arguemnts
     p_order = optval(order, 10)

     n = size(A,1)
     allocate(A2(1:n,1:n)); allocate(X(1:n,1:n))
     allocate(Q(1:n,1:n));  allocate(invQ(1:n,1:n))

     !> Compute the L-infinity norm.
     a_norm = norm_linf(A)
     !> Determine a scaling factor for the matrix.
     ee = int(log2(a_norm)) + 1
     s  = max(0, ee + 1)
     !> Scale the input matrix & initialize polynomial
     A2 = A / 2.0_wp**s
     X = A2

     !> Initialize P & Q and add first step
     c = 0.5_wp
     E = 0.0_wp; forall (k=1:n) E(k, k) = 1.0_wp
     E = E + c*A2
   
     Q = 0.0_wp; forall (k=1:n) Q(k, k) = 1.0_wp
     Q = Q - c*A2
   
     !> Iteratively compute Pade approximation n,n
     p = .true.
     do k = 2, p_order
       c = c*( p_order - k + 1 )/( k*( 2*p_order - k + 1 ) )
       X = matmul( A2, X )
       E = E + c*X
       if ( p ) then
         Q = Q + c*X
       else
         Q = Q - c*X
       end if
       p = .not. p
     end do

     !>  E -> Q^(-1) * P
     invQ = Q
     call inv(invQ)
     E = matmul(invQ,E)

     !>  E -> E^(2*S)
     do k = 1, s
       E = matmul ( E, E )
     end do
   
     return
   end subroutine expm

   function log2(x) result(y)
     ! compute the base-2 logarithm of the input
     implicit none
     real(kind=wp), intent(in) :: x
     real(kind=wp) :: y
     y = log(x) / log(2.0_wp)
   end function log2

   function norm_linf(A) result(norm)
      ! compute the infinity norm of the real-valued input matrix A
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

   function norm_fro(A) result(norm)
      ! compute the frobenuius norm of the real-valued input matrix A
      implicit none   
      real(kind=wp), intent(in) :: A(:,:)
      !> internals
      integer       :: i, j, n
      real(kind=wp) :: norm
      !> initialize
      norm = 0.0_wp
      n = size(A,1)
      do i = 1, n
         do j = 1, n
            norm = norm + A(i,j)**2
         end do
      end do
      norm = sqrt(norm)
   end function norm_fro

   function factorial(n) result(y)
      ! compute n!
      implicit none
      integer, intent(in) :: n 
      integer :: i, y
      y = 1
      do i = 2, n
         y = y*i 
      end do
   end function factorial

end module LightROM_expmlib