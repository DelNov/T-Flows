!==============================================================================!
  subroutine Wall_Time_Max_Hours(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads the maximum number of wall-clock hours of the computation from the
!>  control file.  It is useful when you run your simulation in queues with
!>  time limits.  Default is 168, which is one week.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real                :: val      !! maximim number od wall-clock hours
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  ! 168 hours is one week
  call Control % Read_Real_Item('WALL_TIME_MAX_HOURS', 168.0, val, verbose)

  end subroutine
