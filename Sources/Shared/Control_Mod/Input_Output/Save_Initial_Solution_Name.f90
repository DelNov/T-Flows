!==============================================================================!
  subroutine Save_Initial_Solution_Name(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  character(SL)       :: val
  logical,   optional :: verbose
!==============================================================================!

  call Control % Read_Char_Item('SAVE_INITIAL_SOLUTION_NAME', 'skip',  &
                                 val, verbose)

  end subroutine
