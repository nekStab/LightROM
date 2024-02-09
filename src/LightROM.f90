module LightROM
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, LightROM!"
  end subroutine say_hello
end module LightROM
