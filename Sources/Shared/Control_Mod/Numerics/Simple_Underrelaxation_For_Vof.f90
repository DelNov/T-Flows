!==============================================================================!
  subroutine Control_Mod_Simple_Underrelaxation_For_Vof(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('SIMPLE_UNDERRELAXATION_FOR_VOF',  &
                                   0.5, val, verbose)

  end subroutine
