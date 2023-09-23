!==============================================================================!
  pure subroutine Set_First_Dt(Time, val)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Time_Type), intent(inout) :: Time
  integer,          intent(in)    :: val
!==============================================================================!

  Time % first_time_step = val

  end subroutine