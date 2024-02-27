program Tester
  !> Fortran best practices.
  use, intrinsic :: iso_fortran_env, only : error_unit
  !> Unit-test utility.
  use testdrive, only : run_testsuite, new_testsuite, testsuite_type
  !> Only dummy test. Needs to be removed later.
  use testdrive, only : new_unittest, unittest_type, error_type, check
  !> Abstract implementation of ROM-LTI techniques.
  use LightROM
  use TestVector
  use TestMatrices
  use TestExpm
  use TestLyapunov

  implicit none

  !> Unit-test related.
  integer :: status, i
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("+", *(1x, a))'

  !-------------------------------
  !-----                     -----
  !-----     BEGIN TESTS     -----
  !-----                     -----
  !-------------------------------

  !> Display information about the version of LightKrylov being tested.
  call greetings()

  !> Test status.
  status = 0

  !> List of test suites.
  testsuites = [new_testsuite("Matrix Exponential", collect_expm_testsuite)]
  !testsuites = [new_testsuite("Lyapunov Utils", collect_lyapunov_utils_testsuite)]

  !> Run each test suite.
  do i = 1, size(testsuites)
     write(*, *)            "------------------------------"
     write(error_unit, fmt) "Testing :", testsuites(i)%name
     write(*, *)            "------------------------------"
     write(*, *)
     call run_testsuite(testsuites(i)%collect, error_unit, status)
     write(*, *)
  enddo

  !> Wrap-up.
  if (status > 0) then
     write(error_unit, '(i0, 1x, a)') status, "test(s) failed!"
     error stop
  else if (status == 0) then
     write(*, *) "All test successfully passed!"
  endif

contains

  subroutine collect_dummy_testsuite(testsuite)
    !> Collection of tests.
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [new_unittest("Dummy test 1", dummy_test_1)]

    return
  end subroutine collect_dummy_testsuite

  subroutine dummy_test_1(error)
    !> Error-type to be returned.
    type(error_type), allocatable, intent(out) :: error

    !> Check if 1 == 1.
    call check(error, 1 == 1)
    return
  end subroutine dummy_test_1

end program Tester
