!==============================================================================!
  subroutine Wall_Time_Max_Hours(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real                :: val
  logical,   optional :: verbose
!==============================================================================!

  ! 168 hours is one week
  call Control % Read_Real_Item('WALL_TIME_MAX_HOURS', 168.0, val, verbose)

  end subroutine
