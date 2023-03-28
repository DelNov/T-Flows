!==============================================================================!
  subroutine Update_By_Rank(Prof, i_fun)
!------------------------------------------------------------------------------!
!   Updates function's timer by her rank (number)                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Profiler_Type), target :: Prof
  integer, intent(in)          :: i_fun
!-----------------------------------[Locals]-----------------------------------!
  real(DP) :: wall_time
!==============================================================================!

  ! Store the last time which was recorded
  Prof % i_time_prev = Prof % i_time_curr

  ! Refresh the value of time_curr
  call system_clock(Prof % i_time_curr)

  ! Calculate how much wall time has passed
  wall_time = real(Prof % i_time_curr - Prof % i_time_prev)  &
            / real(Prof % sys_count_rate)

  if(i_fun > 0) then
    Prof % funct_time(i_fun) = Prof % funct_time(i_fun) + wall_time
  end if

  end subroutine
