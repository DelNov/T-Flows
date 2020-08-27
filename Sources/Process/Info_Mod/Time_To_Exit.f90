!==============================================================================!
  logical function Info_Mod_Time_To_Exit()
!------------------------------------------------------------------------------!
!  Returns true if it is time to exit Process because maximum wall clock time  !
!  has been reached                                                            !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  Info_Mod_Time_To_Exit = real(sys_clock % cur - sys_clock % ini)  &
                        / real(sys_clock % cnt)                    &
                        > sys_clock % wall_time_max

  end function
