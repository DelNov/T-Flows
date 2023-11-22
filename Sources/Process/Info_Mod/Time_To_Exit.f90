!==============================================================================!
  logical function Time_To_Exit(Info)
!------------------------------------------------------------------------------!
!  Returns true if it is time to exit Process because maximum wall clock time  !
!  has been reached                                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Info_Type) :: Info
!==============================================================================!

  Time_To_Exit = real(Info % clock % cur - Info % clock % ini)  &
               / real(Info % clock % cnt)                       &
               > Info % clock % wall_time_max

  end function
