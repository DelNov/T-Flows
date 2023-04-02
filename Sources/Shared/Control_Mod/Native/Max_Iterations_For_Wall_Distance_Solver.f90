!==============================================================================!
  subroutine Control_Mod_Max_Iterations_For_Wall_Distance_Solver(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control % Read_Int_Item('MAX_ITERATIONS_FOR_WALL_DISTANCE_SOLVER',  &
                                120, val, verbose)

  end subroutine
