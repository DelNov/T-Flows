!==============================================================================!
  subroutine Control_Mod_Species_Diffusivity(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('SPECIES_DIFFUSIVITY', 1.0e-6, val, verbose)

  end subroutine
