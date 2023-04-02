!==============================================================================!
  subroutine Tolerance_For_Wall_Distance_Solver(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real,   intent(out) :: val
  logical,   optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('TOLERANCE_FOR_WALL_DISTANCE_SOLVER',  &
                                 1.0e-6, val, verbose)

  end subroutine
