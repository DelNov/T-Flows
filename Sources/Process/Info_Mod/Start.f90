!==============================================================================!
  subroutine Info_Mod_Start()
!------------------------------------------------------------------------------!
!  Start Info_Mod by taking system count rate and initial count rate           !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Get system clock clock rate and initial clock count
  call system_clock(count_rate = sys_clock % cnt)
  call system_clock(sys_clock % ini)

  ! Read maximum wall clock hours
  call Control_Mod_Wall_Time_Max_Hours(sys_clock % wall_time_max,  &
                                       verbose=.true.)

  ! Make it in seconds
  sys_clock % wall_time_max = sys_clock % wall_time_max * 3600

  end subroutine
