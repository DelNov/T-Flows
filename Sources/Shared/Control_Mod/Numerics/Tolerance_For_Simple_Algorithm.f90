!==============================================================================!
  subroutine Control_Mod_Tolerance_For_Simple_Algorithm(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(out) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('TOLERANCE_FOR_SIMPLE_ALGORITHM',  &
                                 1.0e-4, val, verbose)

  end subroutine
