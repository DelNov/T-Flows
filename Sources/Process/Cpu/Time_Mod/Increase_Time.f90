!==============================================================================!
  pure subroutine Increase_Time(Time, val)
!------------------------------------------------------------------------------!
!>  Updates the physical time based on the current time step.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Time_Type), intent(inout) :: Time  !! parent, singleton object Time
  real,             intent(in)    :: val   !! physical time increase, addition
                                           !! of the next time step
!==============================================================================!

  Time % physical_time = Time % physical_time + val

  end subroutine
