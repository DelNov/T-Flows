!==============================================================================!
  pure subroutine Set_Time(Time, val)
!---------------------------------[Arguments]----------------------------------!
  class(Time_Type), intent(inout) :: Time
  real,             intent(in)    :: val
!==============================================================================!

  Time % physical_time = val

  end subroutine
