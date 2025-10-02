!==============================================================================!
  logical function Time_To_Exit(Info)
!------------------------------------------------------------------------------!
!>  Checks the current wall-clock time against the maximum allowed time. If the
!>  elapsed wall-clock time exceeds the maximum, the function returns true,
!>  indicating that the simulation should terminate.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Info_Type) :: Info  !! parent, singleton object Info
!==============================================================================!

  Time_To_Exit = real(Info % clock % cur - Info % clock % ini)  &
               / real(Info % clock % cnt)                       &
               > Info % clock % wall_time_max

  end function
