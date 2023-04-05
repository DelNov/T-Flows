!==============================================================================!
  subroutine Max_Smoothing_Cycles_Curvature_Vof(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control % Read_Int_Item('MAX_SMOOTHING_CYCLES_CURVATURE_VOF',  &
                                2, val, verbose)

  end subroutine
