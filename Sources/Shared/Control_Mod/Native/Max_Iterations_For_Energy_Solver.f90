!==============================================================================!
  subroutine Control_Mod_Max_Iterations_For_Energy_Solver(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control % Read_Int_Item('MAX_ITERATIONS_FOR_ENERGY_SOLVER',  &
                                5, val, verbose)

  end subroutine
