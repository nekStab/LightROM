module lightrom_AbstractLTIsystems
  !> Use the abstract linear operator types defined in LightKrylov.
  use LightKrylov, only : abstract_linop
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
     !> Dynamics matrix.
     class(abstract_linop), allocatable :: A
     !> Input-to-state matrix.
     class(abstract_linop), allocatable :: B
     !> State-to-output matrix.
     class(abstract_linop), allocatable :: C
     !> Feedthrough matrix.
     real(kind=wp)        , allocatable :: D(:, :)
     contains
       private
  end type abstract_lti_system

  !> Abstract discrete LTI system.
  type, extends(abstract_dynamical_system), abstract, public :: abstract_dlti_system
     !> Dynamic matrix.
     class(abstract_linop), allocatable :: A
     !> Input-to-state matrix.
     class(abstract_linop), allocatable :: B
     !> State-to-output matrix.
     class(abstract_linop), allocatable :: C
     !> Feedthrough matrix.
     real(kind=wp)        , allocatable :: D(:, :)
     !> Sampling period.
     real(kind=wp)                      :: dt = 1.0_wp
   contains
     private
  end type abstract_dlti_system

contains

end module lightrom_AbstractLTIsystems
