!==============================================================================!
  subroutine Start_Info(Info)
!------------------------------------------------------------------------------!
!>  Initializes the Info_Mod module at the start of the simulation. It retrieves
!>  the system clock's count rate and the initial count to accurately track the
!>  elapsed time during the simulation. Additionally, it reads and sets the
!>  maximum allowable wall-clock time for the simulation from the control file,
!>  converting this time into seconds for consistent time tracking.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Info_Type) :: Info  !! parent, singleton object Info
!==============================================================================!

  ! Get system clock clock rate and initial clock count
  call system_clock(count_rate = Info % clock % cnt)
  call system_clock(Info % clock % ini)

  ! Read maximum wall clock hours
  call Control % Wall_Time_Max_Hours(Info % clock % wall_time_max,  &
                                     verbose=.true.)

  ! Make it in seconds
  Info % clock % wall_time_max = Info % clock % wall_time_max * 3600

  end subroutine
