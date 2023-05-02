!==============================================================================!
  subroutine Hybrid_Les_Rans_Switch(Control, val, verbose)
!------------------------------------------------------------------------------!
!   Reading switch for hybrid LES/RANS model.                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)        :: Control
  character(SL), intent(out) :: val
  logical,       optional    :: verbose
!==============================================================================!

  call Control % Read_Char_Item('HYBRID_LES_RANS_SWITCH',   &
                                'SWITCH_DISTANCE',          &
                                 val, verbose)
  call String % To_Upper_Case(val)

  end subroutine
