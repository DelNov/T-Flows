!==============================================================================!
  pure integer function Curr_Dt(Time)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Time_Type), intent(in) :: Time
!==============================================================================!

  Curr_Dt = Time % current_time_step

  end function
