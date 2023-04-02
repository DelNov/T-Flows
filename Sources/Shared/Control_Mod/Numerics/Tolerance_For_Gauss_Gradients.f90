!==============================================================================!
  subroutine Control_Mod_Tolerance_For_Gauss_Gradients(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(out) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('TOLERANCE_FOR_GAUSS_GRADIENTS',  &
                                 1.0e-3, val, verbose)

  end subroutine
