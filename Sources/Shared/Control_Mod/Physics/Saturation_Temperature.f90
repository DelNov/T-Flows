!==============================================================================!
  subroutine Control_Mod_Saturation_Temperature(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(out) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('SATURATION_TEMPERATURE', 100., val, verbose)

  end subroutine
