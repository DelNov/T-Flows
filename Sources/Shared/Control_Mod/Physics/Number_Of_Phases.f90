!==============================================================================!
  subroutine Number_Of_Phases(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control % Read_Int_Item('NUMBER_OF_PHASES', 1, val, verbose)

  end subroutine
