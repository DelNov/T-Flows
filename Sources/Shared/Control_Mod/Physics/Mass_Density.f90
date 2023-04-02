!==============================================================================!
  subroutine Control_Mod_Mass_Density(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(out) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('MASS_DENSITY', 1.0, val, verbose)

  end subroutine
