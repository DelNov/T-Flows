!==============================================================================!
  subroutine Control_Mod_Tolerance_For_Turbulence_Solver(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(out) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('TOLERANCE_FOR_TURBULENCE_SOLVER',  &
                                 1.0e-6, val, verbose)

  end subroutine
