!==============================================================================!
  subroutine Control_Mod_Latent_Heat(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(out) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('LATENT_HEAT', 1.0, val, verbose)

  end subroutine
