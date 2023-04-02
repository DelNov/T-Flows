!==============================================================================!
  subroutine Surface_Tension(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real,   intent(out) :: val
  logical,   optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('SURFACE_TENSION', 0.0, val, verbose)

  end subroutine
