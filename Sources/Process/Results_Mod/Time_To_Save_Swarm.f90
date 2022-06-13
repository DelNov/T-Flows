!==============================================================================!
  logical function Time_To_Save_Swarm(Results, curr_dt)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Results_Type) :: Results
  integer             :: curr_dt  ! current time step
!==============================================================================!

  Time_To_Save_Swarm = mod(curr_dt, Results % interval_swarm) .eq. 0  &
                       .or.                                           &
                       curr_dt .eq. 0 .and. Results % initial

  end function
