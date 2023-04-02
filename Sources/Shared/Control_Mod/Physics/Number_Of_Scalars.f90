!==============================================================================!
  subroutine Number_Of_Scalars(Control, val, verbose)
!------------------------------------------------------------------------------!
!   Reading stuff related to passive scalars                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control % Read_Int_Item('NUMBER_OF_SCALARS', 0, val, verbose)

  end subroutine
