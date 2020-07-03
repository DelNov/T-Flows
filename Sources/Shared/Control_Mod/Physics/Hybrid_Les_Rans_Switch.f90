!==============================================================================!
  subroutine Control_Mod_Hybrid_Les_Rans_Switch(val, verbose)
!------------------------------------------------------------------------------!
!   Reading switch for hybrid LES/RANS model.                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=80) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('HYBRID_LES_RANS_SWITCH',   &
                                  'SWITCH_DISTANCE',          &
                                   val, verbose)
  call To_Upper_Case(val)

  end subroutine
