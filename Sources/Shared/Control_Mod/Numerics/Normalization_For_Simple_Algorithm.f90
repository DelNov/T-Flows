!==============================================================================!
  subroutine Control_Mod_Normalization_For_Simple_Algorithm(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('NORMALIZATION_FOR_SIMPLE_ALGORITHM',  &
                                   1.0, val, verbose)

  end subroutine
