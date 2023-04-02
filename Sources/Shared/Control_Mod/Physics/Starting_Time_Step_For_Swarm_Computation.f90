!==============================================================================!
  subroutine Starting_Time_Step_For_Swarm_Computation(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control % Read_Int_Item('STARTING_TIME_STEP_FOR_SWARM_COMPUTATION',  &
                                1200, val, verbose)

  end subroutine
