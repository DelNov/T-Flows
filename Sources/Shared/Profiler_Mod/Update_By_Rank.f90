!==============================================================================!
  subroutine Update_By_Rank(Prof, i_fun)
!------------------------------------------------------------------------------!
!>  Updates the accumulated time for a specified function. It calculates the
!>  elapsed time since the last update and adds this to the total time spent
!>  in the function.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Profiler_Type), target :: Prof   !! parent class
  integer, intent(in)          :: i_fun  !! function rank
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
