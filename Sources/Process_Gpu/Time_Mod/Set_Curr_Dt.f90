!==============================================================================!
  pure subroutine Set_Curr_Dt(Time, val)
!------------------------------------------------------------------------------!
!>  Function to set the current time step stored in Time.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Time_Type), intent(inout) :: Time  !! parent, singleton object Time
  integer,          intent(in)    :: val   !! desired value of the time step
!==============================================================================!

  Time % current_time_step = val

  end subroutine
