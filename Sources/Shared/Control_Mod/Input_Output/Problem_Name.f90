!==============================================================================!
  subroutine Read_Problem_Name(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  character(SL)       :: val
  logical,   optional :: verbose
!==============================================================================!

  call Control % Read_Char_Item('PROBLEM_NAME', 'unknown', val, verbose)

  end subroutine
