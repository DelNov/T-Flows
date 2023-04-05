!==============================================================================!
  subroutine Tolerance_For_Simple_Algorithm(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real,   intent(out) :: val
  logical,   optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('TOLERANCE_FOR_SIMPLE_ALGORITHM',  &
                                 1.0e-4, val, verbose)

  end subroutine
