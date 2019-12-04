!==============================================================================!
  subroutine Control_Mod_Problem_Name(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=80) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('PROBLEM_NAME', 'unknown',  &
                                   val, verbose)

  end subroutine
