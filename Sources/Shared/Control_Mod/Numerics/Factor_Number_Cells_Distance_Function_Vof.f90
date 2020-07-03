!==============================================================================!
  subroutine Control_Mod_Factor_Number_Cells_Distance_Function_Vof(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('FACTOR_NUMBER_CELLS_DISTANCE_FUNCTION',  &
                                   1.5, val, verbose)

  end subroutine
