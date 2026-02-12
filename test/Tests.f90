program Tester
  !> Fortran best practices.
  use, intrinsic :: iso_fortran_env, only : error_unit
  !> Fortran standard library
  use stdlib_logger, only: information_level, warning_level, debug_level, error_level, none_level
  !> Unit-test utility.
  use testdrive, only : run_testsuite, new_testsuite, testsuite_type
  !> Only dummy test. Needs to be removed later.
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use LightKrylov
  use LightKrylov_Logger
  use LightKrylov_TestUtils
  !> Abstract implementation of ROM-LTI techniques.
  use LightROM
  use LightROM_TestUtils
  use LightROM_TestLyapunov

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
  call greetings_LightROM()

  ! Turn off logging during tests (unless you REALLY want it)
  !call logger%configure(level=error_level, time_stamp=.false.); write(*,*) 'Logging set to error_level.'
  write(*,*) ""; write(*,*) ""

  !> Test status.
  status = 0

  !> List of test suites.
  testsuites = [ & 
               new_testsuite("Lyapunov Utils", collect_lyapunov_utils_testsuite) &
               ]

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

end program Tester
