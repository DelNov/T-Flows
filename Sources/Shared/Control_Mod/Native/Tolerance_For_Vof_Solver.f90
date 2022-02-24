!==============================================================================!
  subroutine Control_Mod_Tolerance_For_Vof_Solver(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(out) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('TOLERANCE_FOR_VOF_SOLVER',  &
                                   1.0e-6, val, verbose)

  end subroutine
