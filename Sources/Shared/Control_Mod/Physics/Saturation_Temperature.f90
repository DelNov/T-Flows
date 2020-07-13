!==============================================================================!
  subroutine Control_Mod_Saturation_Temperature(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('SATURATION_TEMPERATURE', 1.0, val, verbose)

  end subroutine
