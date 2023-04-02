!==============================================================================!
  subroutine Control_Mod_Tolerance_For_Potential_Solver(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(out) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('TOLERANCE_FOR_POTENTIAL_SOLVER',  &
                                 1.0e-6, val, verbose)

  end subroutine
