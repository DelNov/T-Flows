!==============================================================================!
  subroutine Control_Mod_Dynamic_Viscosity(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('LES_DYNAMIC_VISCOSITY', 0.01, val, verbose)

  end subroutine
