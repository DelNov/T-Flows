!==============================================================================!
  subroutine Max_Correction_Cycles_Beta_Vof(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control % Read_Int_Item('MAX_CORRECTION_CYCLES_BETA_VOF',  &
                                2, val, verbose)

  end subroutine
