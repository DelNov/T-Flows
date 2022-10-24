!==============================================================================!
  subroutine Update_By_Rank(Profiler, i_fun)
!------------------------------------------------------------------------------!
!   Updates function's timer by her rank (number)                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Profiler_Type), target :: Profiler
  integer, intent(in)          :: i_fun
!-----------------------------------[Locals]-----------------------------------!
  real(DP) :: wall_time
!==============================================================================!

  ! Store the last time which was recorded
  Profiler % i_time_prev = Profiler % i_time_curr

  ! Refresh the value of time_curr
  call system_clock(Profiler % i_time_curr)

  ! Calculate how much wall time has passed
  wall_time = real(Profiler % i_time_curr - Profiler % i_time_prev)  &
            / real(Profiler % sys_count_rate)

  if(i_fun > 0) then
    Profiler % funct_time(i_fun) = Profiler % funct_time(i_fun) + wall_time
  end if

  end subroutine
