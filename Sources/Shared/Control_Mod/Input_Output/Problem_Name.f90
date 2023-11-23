!==============================================================================!
  subroutine Read_Problem_Name(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads the problem name from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  character(SL)       :: val      !! problem name (the name of the grid you use)
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Char_Item('PROBLEM_NAME', 'unknown', val, verbose)

  end subroutine
