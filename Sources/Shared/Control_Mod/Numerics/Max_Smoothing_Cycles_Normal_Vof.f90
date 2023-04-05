!==============================================================================!
  subroutine Max_Smoothing_Cycles_Normal_Vof(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control % Read_Int_Item('MAX_SMOOTHING_CYCLES_NORMAL_VOF',  &
                                0, val, verbose)

  end subroutine
