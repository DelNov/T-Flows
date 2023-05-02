!==============================================================================!
  subroutine Saturation_Temperature(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real,   intent(out) :: val
  logical,   optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('SATURATION_TEMPERATURE', 100., val, verbose)

  end subroutine
