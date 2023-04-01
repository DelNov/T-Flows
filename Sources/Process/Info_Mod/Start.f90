!==============================================================================!
  subroutine Info_Mod_Start_Info()
!------------------------------------------------------------------------------!
!  Start Info_Mod by taking system count rate and initial count rate           !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Get system clock clock rate and initial clock count
  call system_clock(count_rate = Info % clock % cnt)
  call system_clock(Info % clock % ini)

  ! Read maximum wall clock hours
  call Control_Mod_Wall_Time_Max_Hours(Info % clock % wall_time_max,  &
                                       verbose=.true.)

  ! Make it in seconds
  Info % clock % wall_time_max = Info % clock % wall_time_max * 3600

  end subroutine
