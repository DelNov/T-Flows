!==============================================================================!
  pure subroutine Set_Time(Time, val)
!------------------------------------------------------------------------------!
!>  Function to set the current physical time.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Time_Type), intent(inout) :: Time  !! parent, singleton object Time
  real,             intent(in)    :: val   !! physical time you want to set
!==============================================================================!

  Time % physical_time = val

  end subroutine
