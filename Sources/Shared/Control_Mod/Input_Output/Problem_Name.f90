!==============================================================================!
  subroutine Control_Mod_Problem_Name(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL)     :: val
  logical, optional :: verbose
!==============================================================================!

  call Control % Read_Char_Item('PROBLEM_NAME', 'unknown', val, verbose)

  end subroutine
