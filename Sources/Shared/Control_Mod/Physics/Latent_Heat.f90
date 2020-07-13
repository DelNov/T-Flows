!==============================================================================!
  subroutine Control_Mod_Latent_Heat(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('LATENT_HEAT', 1.0, val, verbose)

  end subroutine
