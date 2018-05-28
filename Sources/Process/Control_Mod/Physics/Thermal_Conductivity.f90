!==============================================================================!
  subroutine Control_Mod_Thermal_Conductivity(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('THERMAL_CONDUCTIVITY', 1.0, val, verbose)

  end subroutine
