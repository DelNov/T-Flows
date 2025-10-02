!==============================================================================!
  subroutine Save_Initial_Solution_Name(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads the name of the file with initial conditions to save from the control
!>  file. (I believe this is not used at all.)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  character(SL)       :: val      !! initial solution name
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Char_Item('SAVE_INITIAL_SOLUTION_NAME', 'skip',  &
                                 val, verbose)

  end subroutine
