!==============================================================================!
  pure subroutine Set_Curr_Dt(Time, val)
!---------------------------------[Arguments]----------------------------------!
  class(Time_Type), intent(inout) :: Time
  integer,          intent(in)    :: val
!==============================================================================!

  Time % current_time_step = val

  end subroutine
