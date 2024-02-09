module lightrom_AbstractLTIsystems
  use LightKrylov
  implicit none
  include "dtypes.h"

  private

  !-------------------------------------------------------
  !-----     ABSTRACT LTI SYSTEM TYPE DEFINITION     -----
  !-------------------------------------------------------

  !> General abstract type for general system.
  type, abstract, public :: abstract_dynamical_system
  end type abstract_dynamical_system

  !> Abstract continuous LTI system.
  type, extends(abstract_dynamical_system), abstract, public :: abstract_lti_system
  end type abstract_lti_system

  !> Abstract discrete LTI system.
  type, extends(abstract_dynamical_system), abstract, public :: abstract_dlti_system
  end type abstract_dlti_system

contains
end module lightrom_AbstractLTIsystems
