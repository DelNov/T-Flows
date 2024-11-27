!==============================================================================!
  pure real function Get_Time(Time)
!------------------------------------------------------------------------------!
!>  Function to retrieve the current physical time of the simulation.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Time_Type), intent(in) :: Time  !! parent, singleton object Time
!==============================================================================!

  Get_Time = Time % physical_time

  end function
