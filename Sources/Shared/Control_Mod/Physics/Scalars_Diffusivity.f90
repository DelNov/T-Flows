!==============================================================================!
  subroutine Control_Mod_Scalars_Diffusivity(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('SCALARS_DIFFUSIVITY', 1.0e-6, val, verbose)

  end subroutine
