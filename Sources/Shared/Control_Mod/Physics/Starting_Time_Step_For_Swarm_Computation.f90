!==============================================================================!
  subroutine Control_Mod_Starting_Time_Step_For_Swarm_Computation(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control_Mod_Read_Int_Item('STARTING_TIME_STEP_FOR_SWARM_COMPUTATION',  &
                                 1200, val, verbose)

  end subroutine
