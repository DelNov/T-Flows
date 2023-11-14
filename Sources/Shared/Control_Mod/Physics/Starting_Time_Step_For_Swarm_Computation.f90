!==============================================================================!
  subroutine Starting_Time_Step_For_Swarm_Computation(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads starting time step for swarm (particles) simulations.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  integer, intent(out) :: val      !! starting time step for swarm computation
  logical, optional    :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Int_Item('STARTING_TIME_STEP_FOR_SWARM_COMPUTATION',  &
                                1200, val, verbose)

  end subroutine
