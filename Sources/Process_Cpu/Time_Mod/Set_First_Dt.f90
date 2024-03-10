!==============================================================================!
  pure subroutine Set_First_Dt(Time, val)
!------------------------------------------------------------------------------!
!>  Function to set the firts time step stored in Time.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Time_Type), intent(inout) :: Time  !! parent, singleton object Time
  integer,          intent(in)    :: val   !! value for the first time step
!==============================================================================!

  Time % first_time_step = val

  end subroutine
