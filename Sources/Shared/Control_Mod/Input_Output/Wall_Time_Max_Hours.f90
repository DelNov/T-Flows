!==============================================================================!
  subroutine Control_Mod_Wall_Time_Max_Hours(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  ! 168 hours is one week
  call Control % Read_Real_Item('WALL_TIME_MAX_HOURS', 168.0, val, verbose)

  end subroutine
