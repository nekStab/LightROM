module TestMatrices
  use LightKrylov
  use TestVector
  
  implicit none

  private

  !---------------------------------------
  !-----     GENERAL REAL MATRIX     -----
  !---------------------------------------
  type, extends(abstract_linop), public :: rmatrix
     real(kind=wp), dimension(test_size, test_size) :: data = 0.0_wp
   contains
     private
     procedure, pass(self), public :: matvec => general_matvec
     procedure, pass(self), public :: rmatvec => general_rmatvec
  end type rmatrix

  !----------------------------------------
  !-----     SYM. POS. DEF MATRIX     -----
  !----------------------------------------
  type, extends(abstract_spd_linop), public :: spd_matrix
     real(kind=wp), dimension(test_size, test_size) :: data = 0.0_wp
   contains
     private
     procedure, pass(self), public :: matvec => spd_matvec
     procedure, pass(self), public :: rmatvec => spd_matvec
  end type spd_matrix

contains

  !-------------------------------------------------------------------------------------
  !-----                                                                           -----
  !-----     DEFINITION OF THE TYPE-BOUND PROCEDURES FOR GENERAL REAL MATRICES     -----
  !-----                                                                           -----
  !-------------------------------------------------------------------------------------

  subroutine general_matvec(self, vec_in, vec_out)
    class(rmatrix)        , intent(in)  :: self
    class(abstract_vector), intent(in)  :: vec_in
    class(abstract_vector), intent(out) :: vec_out

    select type(vec_in)
    type is(rvector)
       select type(vec_out)
       type is(rvector)
          vec_out%data = matmul(self%data, vec_in%data)
       end select
    end select
  end subroutine general_matvec

  subroutine general_rmatvec(self, vec_in, vec_out)
    class(rmatrix), intent(in)          :: self
    class(abstract_vector), intent(in)  :: vec_in
    class(abstract_vector), intent(out) :: vec_out

    select type(vec_in)
    type is(rvector)
       select type(vec_out)
       type is(rvector)
          vec_out%data = matmul(transpose(self%data), vec_in%data)
       end select
    end select
  end subroutine general_rmatvec

  !-----------------------------------------------------------------------------------------------
  !-----                                                                                     -----
  !-----     DEFINITION OF THE TYPE-BOUND PROCEDURES FOR SYM. POS. DEF. LINEAR OPERATORS     -----
  !-----                                                                                     -----
  !-----------------------------------------------------------------------------------------------

  subroutine spd_matvec(self, vec_in, vec_out)
    class(spd_matrix)     , intent(in)  :: self
    class(abstract_vector), intent(in)  :: vec_in
    class(abstract_vector), intent(out) :: vec_out

    select type(vec_in)
    type is(rvector)
       select type(vec_out)
       type is(rvector)
          vec_out%data = matmul(self%data, vec_in%data)
       end select
    end select
    return
  end subroutine spd_matvec

  !----------------------------------------------------------------------------------
  !-----                                                                        -----
  !-----     DEFINITION OF THE TYPE-BOUND PROCEDURES FOR HERMITIAN MATRICES     -----
  !-----                                                                        -----
  !----------------------------------------------------------------------------------

end module TestMatrices
