!==============================================================================!
  subroutine Thermal_Conductivity(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real,   intent(out) :: val
  logical,   optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('THERMAL_CONDUCTIVITY', 1.0, val, verbose)

  end subroutine
