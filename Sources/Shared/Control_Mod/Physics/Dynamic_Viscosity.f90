!==============================================================================!
  subroutine Control_Mod_Dynamic_Viscosity(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(out) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('DYNAMIC_VISCOSITY', 0.01, val, verbose)

  end subroutine
