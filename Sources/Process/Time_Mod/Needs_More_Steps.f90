!==============================================================================!
  logical function Needs_More_Steps(Time)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Time_Type), intent(inout) :: Time
!-----------------------------------[Locals]-----------------------------------!
  logical, save :: first_entry = .true.
!==============================================================================!

  !--------------------------------------------!
  !   You visit this for the first time, you   !
  !    are at the beginning of a simulation    !
  !--------------------------------------------!
  if(first_entry) then
    Time % current_time_step = Time % First_Time_Step + 1
    first_entry = .false.

  !-----------------------------------------!
  !   You are not in the first time step,   !
  !     increase the time step counter      !
  !-----------------------------------------!
  else
    Time % current_time_step =  &
    Time % current_time_step + 1

  end if

  !-------------------------------------------------!
  !   Decide if more time steps are needed or not   !
  !-------------------------------------------------!
  if(Time % current_time_step .gt. Time % last_time_step) then
    Needs_More_Steps = .false.

  else
    Needs_More_Steps = .true.

  end if

  end function
