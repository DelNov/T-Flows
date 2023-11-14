!==============================================================================!
  subroutine Hybrid_Les_Rans_Switch(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads switch for hybrid LES/RANS model from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)        :: Control  !! parent class
  character(SL), intent(out) :: val      !! hybrid LES/RANS switch
  logical,       optional    :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Char_Item('HYBRID_LES_RANS_SWITCH',   &
                                'SWITCH_DISTANCE',          &
                                 val, verbose)
  call String % To_Upper_Case(val)

  end subroutine
