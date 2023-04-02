!==============================================================================!
  subroutine Time_Step(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real,   intent(out) :: val
  logical,   optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('TIME_STEP', 1.0e-2,      &
                                 val, verbose)

  end subroutine
