!==============================================================================!
  logical function Time_To_Save_Swarm(Results)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Results_Type) :: Results  !! parent class
!==============================================================================!

  Time_To_Save_Swarm = mod(Time % Curr_Dt(),                 &
                           Results % interval_swarm) .eq. 0  &
                       .or.                                  &
                       Time % Curr_Dt() .eq. 0 .and. Results % initial

  end function
