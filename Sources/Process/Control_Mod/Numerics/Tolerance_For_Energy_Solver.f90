!==============================================================================!
  subroutine Control_Mod_Tolerance_For_Energy_Solver(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('TOLERANCE_FOR_ENERGY_SOLVER', 1.0e-6,  &
                                   val, verbose)

  end subroutine
