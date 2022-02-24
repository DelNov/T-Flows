!==============================================================================!
  subroutine Control_Mod_Max_Iterations_For_Potential_Solver(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control_Mod_Read_Int_Item('MAX_ITERATIONS_FOR_POTENTIAL_SOLVER',  &
                                  120, val, verbose)

  end subroutine
