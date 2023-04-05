!==============================================================================!
  subroutine Turbulent_Schmidt_Number(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real,   intent(out) :: val
  logical,   optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('TURBULENT_SCHMIDT_NUMBER', 0.9, val, verbose)

  end subroutine
