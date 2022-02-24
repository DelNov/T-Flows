!==============================================================================!
  subroutine Control_Mod_Tolerance_For_Scalars_Solver(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(out) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('TOLERANCE_FOR_SCALARS_SOLVER',  &
                                   1.0e-6, val, verbose)

  end subroutine
