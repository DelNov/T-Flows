!==============================================================================!
  pure integer function Curr_Dt(Time)
!------------------------------------------------------------------------------!
!>  Function to retreive the current time step stored in Time.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Time_Type), intent(in) :: Time  !! parent, singleton object Time
!==============================================================================!

  Curr_Dt = Time % current_time_step

  end function
