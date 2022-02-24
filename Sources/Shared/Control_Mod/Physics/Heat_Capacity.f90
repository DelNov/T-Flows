!==============================================================================!
  subroutine Control_Mod_Heat_Capacity(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(out) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('HEAT_CAPACITY', 1.0, val, verbose)

  end subroutine
