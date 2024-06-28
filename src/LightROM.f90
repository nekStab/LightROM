module LightROM
   ! LTI system definitions
   use LightROM_AbstractLTIsystems
   ! General utilities
   use LightROM_Utils
   ! Solvers for Lyapunov equations
   use LightROM_LyapunovSolvers
   implicit none
   private

   !> Global variables.
   public :: greetings_LightROM

   ! Abstract types
   public :: abstract_lti_system_rdp
   public :: abstract_sym_low_rank_state_rdp

   ! Utils
   public :: dlra_opts

   ! Solvers
   public :: projector_splitting_DLRA_lyapunov_integrator

contains

  subroutine greetings_LightROM()
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
  end subroutine greetings_LightROM

end module LightROM
