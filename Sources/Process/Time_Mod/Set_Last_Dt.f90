!==============================================================================!
  pure subroutine Set_Last_Dt(Time, val)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Time_Type), intent(inout) :: Time
  integer,          intent(in)    :: val
!==============================================================================!

  Time % last_time_step = val

  end subroutine