!==============================================================================!
  logical function Info_Mod_Time_Is_Up()
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  Info_Mod_Time_Is_Up = real(sys_clock % cur - sys_clock % ini)  &
                      / real(sys_clock % cnt)                    &
                      > sys_clock % wall_time_max

  end function
