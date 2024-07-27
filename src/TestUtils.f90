module LightROM_TestUtils
    ! Fortran Standard Library
    ! LightKrylov
    use LightKrylov
    use LightKrylov_Constants
    use LightKrylov_TestUtils
    ! LightROM
    use LightROM_AbstractLTIsystems

    implicit none

    private

    character(len=128), parameter, private :: this_module = 'LightROM_TestUtils'

    !-----------------------------------------------------
    !-----     TEST SYM LOW-RANK STATE DEFINTION     -----
    !-----------------------------------------------------
    
    type, extends(abstract_sym_low_rank_state_rdp), public :: sym_low_rank_state_rdp
    end type sym_low_rank_state_rdp

contains


end module LightROM_TestUtils
