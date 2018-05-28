!==============================================================================!
  subroutine Control_Mod_Blending_Coefficient_Energy(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('BLENDING_COEFFICIENT_ENERGY', 1.0,  &
                                   val, verbose)

  end subroutine
