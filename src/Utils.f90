module lightrom_utils
  use stdlib_optval, only : optval
  implicit none
  include "dtypes.h"

  private
  public :: print_mat

contains

   subroutine print_mat(m,n,A)

      implicit none

      integer , intent(in)                         :: n
      integer , intent(in)                         :: m
      real (kind=wp)  , intent(in)                 :: a(m,n)

      !> internals
      real (kind=wp) :: amax
      real (kind=wp) :: amin
      integer :: i
      character ( len = 10 ) iform
      integer :: ihi
      integer :: ilo
      logical :: integ
      integer :: j
      integer :: jhi
      integer :: jlo
      integer :: lmax
      integer :: npline
!     
!      Check if all entries are integral.
!     
      integ = .true.

      do i = 1, m
        do j = 1, n

          if ( integ ) then
            if ( a(i,j) /= real ( int ( a(i,j) ),kind=wp) ) then
              integ = .false.
            end if
          end if
       
        end do
      end do
!     
!      Find the maximum and minimum entries.
!     
      amax = maxval ( a(1:m,1:n) )
      amin = minval ( a(1:m,1:n) )
!     
!      Use the information about the maximum size of an entry to
!      compute an intelligent format for use with integer entries.
!     
!      Later, we might also do this for real matrices.
!     
      lmax = int ( log10 ( amax ) )
   
      if ( integ ) then
        npline = 79 / ( lmax + 3 )
        write ( iform, '(''('',i2,''I'',i2,'')'')' ) npline, lmax+3
      else
        npline = 10
        iform = ' '
      end if
!     
!      Print a scalar quantity.
!     
      if ( m == 1 .and. n == 1 ) then
      
        if ( integ ) then
          write ( *, iform ) int ( a(1,1) )
        else
          write ( *, '(2x,g10.2)' ) a(1,1)
        end if
!     
!      Column vector of length M,
!     
      else if ( n == 1 ) then
      
        do ilo = 1, m, npline

          ihi = min ( ilo+npline-1, m )

          if ( integ ) then
            write ( *, iform ) ( int ( a(i,1) ), i = ilo, ihi )
          else
            write ( *, '(2x,10g10.2)' ) a(ilo:ihi,1)
          end if
       
        end do
!     
!      Row vector of length N,
!     
      else if ( m == 1 ) then
      
        do jlo = 1, n, npline

          jhi = min ( jlo+npline-1, n )

          if ( integ ) then
            write ( *, iform ) int ( a(1,jlo:jhi) )
          else
            write ( *, '(2x,10g10.2)' ) a(1,jlo:jhi)
          end if
       
        end do
!     
!      M by N Array
!     
      else
      
        do jlo = 1, n, npline

          jhi = min ( jlo+npline-1, n )

          if ( npline < n ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a,i8,a,i8)' ) 'Matrix columns ', jlo, ' to ', jhi
            write ( *, '(a)' ) ' '
          end if
       
          do i = 1, m

            if ( integ ) then
              write ( *, iform ) int ( a(i,jlo:jhi) )
            else
              write ( *, '(2x,5g14.6)' ) a(i,jlo:jhi)
            end if
         
          end do
        end do
     
      end if
      write(*,*)
   
      return
   end subroutine print_mat

end module lightrom_utils
