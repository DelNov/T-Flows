!==============================================================================!
  subroutine Control_Mod_Roughness_Coefficient(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('ROUGHNESS_COEFFICIENT', 0.0,  &
                                   val, verbose)

  end subroutine
