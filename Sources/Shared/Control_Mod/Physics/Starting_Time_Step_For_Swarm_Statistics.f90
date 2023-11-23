!==============================================================================!
  subroutine Starting_Time_Step_For_Swarm_Statistics(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads starting time step for gathering swarm (particles) statistics.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  integer, intent(out) :: val      !! starting time step for swarm statistics
  logical, optional    :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Int_Item('STARTING_TIME_STEP_FOR_SWARM_STATISTICS',  &
                                HUGE_INT, val, verbose)

  end subroutine
