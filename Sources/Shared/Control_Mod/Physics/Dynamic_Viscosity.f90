!==============================================================================!
  subroutine Dynamic_Viscosity(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real,   intent(out) :: val
  logical,   optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('DYNAMIC_VISCOSITY', 0.01, val, verbose)

  end subroutine
