!==============================================================================!
  logical function Time_To_Save(Results, curr_dt)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Results_Type) :: Results
  integer             :: curr_dt  ! current time step
!==============================================================================!

  Time_To_Save = mod(curr_dt, Results % interval) .eq. 0

  end function
