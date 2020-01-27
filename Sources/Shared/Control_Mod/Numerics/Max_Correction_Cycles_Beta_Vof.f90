!==============================================================================!
  subroutine Control_Mod_Max_Correction_Cycles_Beta_Vof(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Int_Item('MAX_CORRECTION_CYCLES_BETA_VOF',  &
                                   2, val, verbose)

  end subroutine
