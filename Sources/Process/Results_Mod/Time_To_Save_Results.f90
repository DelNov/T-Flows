!==============================================================================!
  logical function Time_To_Save_Results(Results)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Results_Type) :: Results
  integer             :: curr_dt  ! current time step
!==============================================================================!

  Time_To_Save_Results = mod(Time % Curr_Dt(), Results % interval) .eq. 0  &
                             .or.                                    &
                             Time % Curr_Dt() .eq. 0 .and. Results % initial

  end function
