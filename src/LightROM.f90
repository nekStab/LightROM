module LightROM
  implicit none
  include "dtypes.h"

  private

  !> Global variables.
  public :: greetings, wp, atol, rtol

contains

  subroutine greetings()
    write(*, *)
    write(*, *)
    write(*, *) "-------------------------------------------------"
    write(*, *) "-------------------------------------------------"
    write(*, *)

    write(*, *) "      _     _       _     _  ______ ________  ___"
    write(*, *) "     | |   (_)     | |   | | | ___ \  _  |  \/  |"
    write(*, *) "     | |    _  __ _| |__ | |_| |_/ / | | | .  . |"
    write(*, *) "     | |   | |/ _` | '_ \| __|    /| | | | |\/| |"
    write(*, *) "     | |___| | (_| | | | | |_| |\ \\ \_/ / |  | |"
    write(*, *) "     \_____/_|\__, |_| |_|\__\_| \_|\___/\_|  |_/"
    write(*, *) "               __/ |"
    write(*, *) "              |___/"
    
    write(*, *)
    write(*, *) "Developped by: Jean-Christophe Loiseau & Simon Kern,"
    write(*, *) "               Arts & Métiers Institute of Technology, 2024,"
    write(*, *) "               jean-christophe.loiseau@ensam.eu"

    write(*, *) "Version -- 0.1.0"
    write(*, *) "License -- BSD 3-Clause"
    write(*, *)

    write(*, *) "-------------------------------------------------"
    write(*, *) "-------------------------------------------------"
    write(*, *)
    write(*, *)
  end subroutine greetings

end module LightROM
