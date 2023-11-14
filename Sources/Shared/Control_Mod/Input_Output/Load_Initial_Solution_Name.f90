!==============================================================================!
  subroutine Load_Initial_Solution_Name(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads the name of the file with initial solution from the control file.
!>  (I believe this is not used at all.)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  character(SL)       :: val      !! name of the file with initial solution
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Char_Item('LOAD_INITIAL_SOLUTION_NAME', 'skip',  &
                                 val, verbose)

  end subroutine
