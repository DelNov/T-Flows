!==============================================================================!
  pure subroutine Set_Last_Dt(Time, val)
!------------------------------------------------------------------------------!
!>  Function to set the last time step of this simulation.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Time_Type), intent(inout) :: Time  !! parent, singleton object Time
  integer,          intent(in)    :: val   !! value for the last time step
!==============================================================================!

  Time % last_time_step = val

  end subroutine
