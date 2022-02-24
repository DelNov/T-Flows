!==============================================================================!
  subroutine Control_Mod_Surface_Tension(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(out) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('SURFACE_TENSION', 0.0, val, verbose)

  end subroutine
