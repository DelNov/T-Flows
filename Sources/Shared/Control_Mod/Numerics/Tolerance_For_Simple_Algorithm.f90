!==============================================================================!
  subroutine Control_Mod_Tolerance_For_Simple_Algorithm(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('TOLERANCE_FOR_SIMPLE_ALGORITHM',  &
                                   1.0e-4, val, verbose)

  end subroutine
