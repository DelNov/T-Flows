!==============================================================================!
  logical function Time_To_Save_Results(Results, curr_dt)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Results_Type) :: Results
  integer             :: curr_dt  ! current time step
!==============================================================================!

  Time_To_Save_Results = mod(curr_dt, Results % interval) .eq. 0     &
                             .or.                                    &
                             curr_dt .eq. 0 .and. Results % initial

  end function
