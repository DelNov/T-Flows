!==============================================================================!
  subroutine Control_Mod_Simple_Underrelaxation_For_Turbulence(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('SIMPLE_UNDERRELAXATION_FOR_TURBULENCE', 0.7,  &
                                   val, verbose)

  end subroutine
