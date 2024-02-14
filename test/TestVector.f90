module TestVector
  use LightKrylov

  implicit none

  private

  public ::  test_size

  integer, parameter :: test_size = 5

  type, extends(abstract_vector), public :: rvector
     real(kind=wp), dimension(test_size) :: data = 0.0_wp
   contains
     private
     procedure, pass(self), public :: zero
     procedure, pass(self), public :: dot
     procedure, pass(self), public :: scal
     procedure, pass(self), public :: axpby
  end type rvector

contains

  !-----------------------------------------------------------
  !-----                                                 -----
  !-----     DEFINITION OF THE TYPE-BOUND PROCEDURES     -----
  !-----                                                 -----
  !-----------------------------------------------------------

  !--> Zero-out a vector.
  subroutine zero(self)
    class(rvector), intent(inout) :: self
    self%data = 0.0_wp
    return
  end subroutine zero

  double precision function dot(self, vec) result(alpha)
    class(rvector), intent(in)         :: self
    class(abstract_vector), intent(in) :: vec

    select type(vec)
    type is(rvector)
       alpha = dot_product(self%data, vec%data)
    end select
    return
  end function dot

  ! --> In-place scalar multiplication.
  subroutine scal(self, alpha)
    class(rvector), intent(inout) :: self
    real(kind=wp), intent(in) :: alpha
    self%data = self%data * alpha
    return
  end subroutine scal

  ! --> axpby interface
  subroutine axpby(self, alpha, vec, beta)
    class(rvector), intent(inout) :: self
    class(abstract_vector), intent(in) :: vec
    real(kind=wp)         , intent(in) :: alpha, beta

    select type(vec)
    type is(rvector)
       self%data = alpha*self%data + beta*vec%data
    end select
    return
  end subroutine axpby

end module TestVector
