!==============================================================================!
  subroutine Smagorinsky_Constant(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real,   intent(out) :: val
  logical,   optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('SMAGORINSKY_CONSTANT', 0.17, val, verbose)

  end subroutine
