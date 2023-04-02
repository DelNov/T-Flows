!==============================================================================!
  subroutine Control_Mod_Turbulent_Schmidt_Number(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(out) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('TURBULENT_SCHMIDT_NUMBER', 0.9, val, verbose)

  end subroutine
