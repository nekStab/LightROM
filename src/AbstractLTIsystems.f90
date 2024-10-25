module LightROM_AbstractLTIsystems
   !> Use the abstract linear operator types defined in LightKrylov.
   use LightKrylov, only : abstract_linop_rdp, abstract_vector_rdp, wp => dp
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
   type, extends(abstract_dynamical_system), abstract, public :: abstract_lti_system_rdp
      !> Dynamics matrix.
      class(abstract_linop_rdp),  allocatable :: A
      !> Exponential propagator.
      class(abstract_linop_rdp),  allocatable :: prop
      !> Input-to-state matrix.
      class(abstract_vector_rdp), allocatable :: B(:)
      !> State-to-output matrix.
      class(abstract_vector_rdp), allocatable :: CT(:)
      !> Feedthrough matrix.
      real(wp),                   allocatable :: D(:, :)
   contains
   end type abstract_lti_system_rdp

   !> Abstract discrete LTI system.
   type, extends(abstract_dynamical_system), abstract, public :: abstract_dlti_system_rdp
      !> Dynamic matrix.
      class(abstract_linop_rdp),  allocatable :: A
      !> Input-to-state matrix.
      class(abstract_vector_rdp), allocatable :: B(:)
      !> State-to-output matrix.
      class(abstract_vector_rdp), allocatable :: CT(:)
      !> Feedthrough matrix.
      real(wp),                   allocatable :: D(:, :)
      !> Sampling period.
      real(wp)                                :: dt = 1.0_wp
   contains
     private
   end type abstract_dlti_system_rdp

   !--------------------------------------------------------------------
   !-----     ABSTRACT LOW RANK REPRESENTATION TYPE DEFINITION     -----
   !--------------------------------------------------------------------

   !> General abstract type for general system.
   type, abstract, public :: abstract_low_rank_representation
   end type abstract_low_rank_representation

   !> Abstract symmetric low-rank representation.
   type, extends(abstract_low_rank_representation), abstract, public :: abstract_sym_low_rank_state_rdp
      !> Low-Rank basis.
      class(abstract_vector_rdp), allocatable :: U(:)
      !> Coefficients
      real(wp),                   allocatable :: S(:, :)
      !> Current approximation rank
      integer                                 :: rk = 1
      !> Has rank been initialized? (for rank-adaptive DLRA)
      logical                                 :: rank_is_initialized = .false.
      !> Simulation time
      real(wp)                                :: time = 0.0_wp
   contains
   end type abstract_sym_low_rank_state_rdp

contains

end module LightROM_AbstractLTIsystems
