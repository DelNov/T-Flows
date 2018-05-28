!==============================================================================!
  subroutine Control_Mod_Turbulent_Prandtl_Number(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('TURBULENT_PRANDTL_NUMBER', 0.9,  &
                                   val, verbose)

  end subroutine
