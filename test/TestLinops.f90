module TestLinops
   ! Fortran Standard Library
   use stdlib_math, only: is_close, all_close

   ! LightKrylov
   use LightKrylov
   use LightKrylov_Constants
   use LightKrylov_Logger
   
   ! Testdrive
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use TestVectors

   implicit none
   
   private

   character*128, parameter, private :: this_module = 'LightKrylov_TestLinops'

   type, extends(abstract_linop_rsp), public :: linop_rsp
       real(sp), dimension(test_size, test_size) :: data = zero_rsp
   contains
       private
       procedure, pass(self), public :: matvec  => matvec_rsp
       procedure, pass(self), public :: rmatvec => rmatvec_rsp
   end type

   type, extends(abstract_sym_linop_rsp), public :: spd_linop_rsp
       real(sp), dimension(test_size, test_size) :: data = zero_rsp
   contains
       private
       procedure, pass(self), public :: matvec => sdp_matvec_rsp
       procedure, pass(self), public :: rmatvec => sdp_matvec_rsp
   end type

   type, extends(abstract_linop_rdp), public :: linop_rdp
       real(dp), dimension(test_size, test_size) :: data = zero_rdp
   contains
       private
       procedure, pass(self), public :: matvec  => matvec_rdp
       procedure, pass(self), public :: rmatvec => rmatvec_rdp
   end type

   type, extends(abstract_sym_linop_rdp), public :: spd_linop_rdp
       real(dp), dimension(test_size, test_size) :: data = zero_rdp
   contains
       private
       procedure, pass(self), public :: matvec => sdp_matvec_rdp
       procedure, pass(self), public :: rmatvec => sdp_matvec_rdp
   end type

   type, extends(abstract_linop_csp), public :: linop_csp
       complex(sp), dimension(test_size, test_size) :: data = zero_csp
   contains
       private
       procedure, pass(self), public :: matvec  => matvec_csp
       procedure, pass(self), public :: rmatvec => rmatvec_csp
   end type

   type, extends(abstract_hermitian_linop_csp), public :: hermitian_linop_csp
       complex(sp), dimension(test_size, test_size) :: data = zero_csp
   contains
       private
       procedure, pass(self), public :: matvec => hermitian_matvec_csp
       procedure, pass(self), public :: rmatvec => hermitian_matvec_csp
   end type

   type, extends(abstract_linop_cdp), public :: linop_cdp
       complex(dp), dimension(test_size, test_size) :: data = zero_cdp
   contains
       private
       procedure, pass(self), public :: matvec  => matvec_cdp
       procedure, pass(self), public :: rmatvec => rmatvec_cdp
   end type

   type, extends(abstract_hermitian_linop_cdp), public :: hermitian_linop_cdp
       complex(dp), dimension(test_size, test_size) :: data = zero_cdp
   contains
       private
       procedure, pass(self), public :: matvec => hermitian_matvec_cdp
       procedure, pass(self), public :: rmatvec => hermitian_matvec_cdp
   end type

contains

   !--------------------------------------------------------------------
   !-----     DEFINITIONS OF THE VARIOUS TYPE-BOUND PROCEDURES     -----
   !--------------------------------------------------------------------

   subroutine matvec_rsp(self, vec_in, vec_out)
       class(linop_rsp), intent(in)  :: self
       class(abstract_vector_rsp)       , intent(in)  :: vec_in
       class(abstract_vector_rsp)       , intent(out) :: vec_out

       select type(vec_in)
       type is(vector_rsp)
           select type(vec_out)
           type is(vector_rsp)

           vec_out%data = matmul(self%data, vec_in%data)

           end select
       end select

       return
   end subroutine matvec_rsp

   subroutine rmatvec_rsp(self, vec_in, vec_out)
       class(linop_rsp), intent(in)  :: self
       class(abstract_vector_rsp)       , intent(in)  :: vec_in
       class(abstract_vector_rsp)       , intent(out) :: vec_out

       select type(vec_in)
       type is(vector_rsp)
           select type(vec_out)
           type is(vector_rsp)

           vec_out%data = matmul(transpose(self%data), vec_in%data)

           end select
       end select

       return
   end subroutine rmatvec_rsp

   subroutine sdp_matvec_rsp(self, vec_in, vec_out)
       class(spd_linop_rsp), intent(in)  :: self
       class(abstract_vector_rsp)       , intent(in)  :: vec_in
       class(abstract_vector_rsp)       , intent(out) :: vec_out

       select type(vec_in)
       type is(vector_rsp)
           select type(vec_out)
           type is(vector_rsp)

           vec_out%data = matmul(self%data, vec_in%data)

           end select
       end select

       return
   end subroutine sdp_matvec_rsp

   subroutine matvec_rdp(self, vec_in, vec_out)
       class(linop_rdp), intent(in)  :: self
       class(abstract_vector_rdp)       , intent(in)  :: vec_in
       class(abstract_vector_rdp)       , intent(out) :: vec_out

       select type(vec_in)
       type is(vector_rdp)
           select type(vec_out)
           type is(vector_rdp)

           vec_out%data = matmul(self%data, vec_in%data)

           end select
       end select

       return
   end subroutine matvec_rdp

   subroutine rmatvec_rdp(self, vec_in, vec_out)
       class(linop_rdp), intent(in)  :: self
       class(abstract_vector_rdp)       , intent(in)  :: vec_in
       class(abstract_vector_rdp)       , intent(out) :: vec_out

       select type(vec_in)
       type is(vector_rdp)
           select type(vec_out)
           type is(vector_rdp)

           vec_out%data = matmul(transpose(self%data), vec_in%data)

           end select
       end select

       return
   end subroutine rmatvec_rdp

   subroutine sdp_matvec_rdp(self, vec_in, vec_out)
       class(spd_linop_rdp), intent(in)  :: self
       class(abstract_vector_rdp)       , intent(in)  :: vec_in
       class(abstract_vector_rdp)       , intent(out) :: vec_out

       select type(vec_in)
       type is(vector_rdp)
           select type(vec_out)
           type is(vector_rdp)

           vec_out%data = matmul(self%data, vec_in%data)

           end select
       end select

       return
   end subroutine sdp_matvec_rdp

   subroutine matvec_csp(self, vec_in, vec_out)
       class(linop_csp), intent(in)  :: self
       class(abstract_vector_csp)       , intent(in)  :: vec_in
       class(abstract_vector_csp)       , intent(out) :: vec_out

       select type(vec_in)
       type is(vector_csp)
           select type(vec_out)
           type is(vector_csp)

           vec_out%data = matmul(self%data, vec_in%data)

           end select
       end select

       return
   end subroutine matvec_csp

   subroutine rmatvec_csp(self, vec_in, vec_out)
       class(linop_csp), intent(in)  :: self
       class(abstract_vector_csp)       , intent(in)  :: vec_in
       class(abstract_vector_csp)       , intent(out) :: vec_out

       select type(vec_in)
       type is(vector_csp)
           select type(vec_out)
           type is(vector_csp)

           vec_out%data = matmul(transpose(conjg(self%data)), vec_in%data)

           end select
       end select

       return
   end subroutine rmatvec_csp

   subroutine hermitian_matvec_csp(self, vec_in, vec_out)
       class(hermitian_linop_csp), intent(in)  :: self
       class(abstract_vector_csp)       , intent(in)  :: vec_in
       class(abstract_vector_csp)       , intent(out) :: vec_out

       select type(vec_in)
       type is(vector_csp)
           select type(vec_out)
           type is(vector_csp)

           vec_out%data = matmul(self%data, vec_in%data)

           end select
       end select

       return
   end subroutine hermitian_matvec_csp

   subroutine matvec_cdp(self, vec_in, vec_out)
       class(linop_cdp), intent(in)  :: self
       class(abstract_vector_cdp)       , intent(in)  :: vec_in
       class(abstract_vector_cdp)       , intent(out) :: vec_out

       select type(vec_in)
       type is(vector_cdp)
           select type(vec_out)
           type is(vector_cdp)

           vec_out%data = matmul(self%data, vec_in%data)

           end select
       end select

       return
   end subroutine matvec_cdp

   subroutine rmatvec_cdp(self, vec_in, vec_out)
       class(linop_cdp), intent(in)  :: self
       class(abstract_vector_cdp)       , intent(in)  :: vec_in
       class(abstract_vector_cdp)       , intent(out) :: vec_out

       select type(vec_in)
       type is(vector_cdp)
           select type(vec_out)
           type is(vector_cdp)

           vec_out%data = matmul(transpose(conjg(self%data)), vec_in%data)

           end select
       end select

       return
   end subroutine rmatvec_cdp

   subroutine hermitian_matvec_cdp(self, vec_in, vec_out)
       class(hermitian_linop_cdp), intent(in)  :: self
       class(abstract_vector_cdp)       , intent(in)  :: vec_in
       class(abstract_vector_cdp)       , intent(out) :: vec_out

       select type(vec_in)
       type is(vector_cdp)
           select type(vec_out)
           type is(vector_cdp)

           vec_out%data = matmul(self%data, vec_in%data)

           end select
       end select

       return
   end subroutine hermitian_matvec_cdp

end module TestLinops

