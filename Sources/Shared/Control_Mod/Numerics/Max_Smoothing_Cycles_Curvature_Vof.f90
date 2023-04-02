!==============================================================================!
  subroutine Control_Mod_Max_Smoothing_Cycles_Curvature_Vof(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control % Read_Int_Item('MAX_SMOOTHING_CYCLES_CURVATURE_VOF',  &
                                2, val, verbose)

  end subroutine
