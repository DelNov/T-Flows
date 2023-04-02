!==============================================================================!
  subroutine Number_Of_Domains(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control % Read_Int_Item('NUMBER_OF_DOMAINS', 1, val, verbose)

  end subroutine
