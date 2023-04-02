!==============================================================================!
  subroutine Control_Mod_Normalization_For_Turbulence_Solver(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('NORMALIZATION_FOR_TURBULENCE_SOLVER',  &
                                 1.0e-6, val, verbose)

  end subroutine
