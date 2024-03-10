!==============================================================================!
  pure integer function First_Dt(Time)
!------------------------------------------------------------------------------!
!>  Function to retreive the first time step of this simulation.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Time_Type), intent(in) :: Time  !! parent, singleton object Time
!==============================================================================!

  First_Dt = Time % first_time_step

  end function
