!==============================================================================!
  pure integer function Last_Dt(Time)
!------------------------------------------------------------------------------!
!>  Function to retreive the last time step of this simulation.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Time_Type), intent(in) :: Time  !! parent, singleton object Time
!==============================================================================!

  Last_Dt = Time % last_time_step

  end function
