!==============================================================================!
  subroutine Control_Mod_Max_Smoothing_Cycles_Normal_Vof(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control % Read_Int_Item('MAX_SMOOTHING_CYCLES_NORMAL_VOF',  &
                                0, val, verbose)

  end subroutine
