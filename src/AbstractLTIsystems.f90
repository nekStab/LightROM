module LightROM_AbstractLTIsystems
   !> Use the abstract linear operator types defined in LightKrylov.
   use LightKrylov, only : abstract_linop, abstract_vector
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
      class(abstract_linop),  allocatable :: A
      !> Input-to-state matrix.
      class(abstract_vector), allocatable :: B(:)
      !> State-to-output matrix.
      class(abstract_vector), allocatable :: CT(:)
      !> Feedthrough matrix.
      real(kind=wp)        ,  allocatable :: D(:, :)
   contains
   end type abstract_lti_system

   !> Abstract discrete LTI system.
   type, extends(abstract_dynamical_system), abstract, public :: abstract_dlti_system
      !> Dynamic matrix.
      class(abstract_linop),  allocatable :: A
      !> Input-to-state matrix.
      class(abstract_vector), allocatable :: B(:)
      !> State-to-output matrix.
      class(abstract_vector), allocatable :: CT(:)
      !> Feedthrough matrix.
      real(kind=wp)        ,  allocatable :: D(:, :)
      !> Sampling period.
      real(kind=wp)                       :: dt = 1.0_wp
   contains
     private
   end type abstract_dlti_system

   !> Abstract reduced-order continuous LTI system.
   type, extends(abstract_dynamical_system), abstract, public :: abstract_ROM_lti_system
      !> Dynamics matrix.
      real(kind=wp)        ,  allocatable :: A(:, :)
      !> Input-to-state matrix.
      real(kind=wp)        ,  allocatable :: B(:, :)
      !> State-to-output matrix.
      real(kind=wp)        ,  allocatable :: C(:, :)
      !> Feedthrough matrix.
      real(kind=wp)        ,  allocatable :: D(:, :)
   contains
   end type abstract_ROM_lti_system

   !--------------------------------------------------------------------
   !-----     ABSTRACT LOW RANK REPRESENTATION TYPE DEFINITION     -----
   !--------------------------------------------------------------------

   !> General abstract type for general system.
   type, abstract, public :: abstract_low_rank_representation
   end type abstract_low_rank_representation

   !> Abstract symmetric low-rank representation.
   type, extends(abstract_low_rank_representation), abstract, public :: abstract_sym_low_rank_state
      !> Low-Rank basis.
      class(abstract_vector), allocatable :: U(:)
      !> Coefficients
      real(kind=wp)        ,  allocatable :: S(:, :)
   contains
   end type abstract_sym_low_rank_state

contains

end module LightROM_AbstractLTIsystems
