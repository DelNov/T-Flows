!==============================================================================!
  subroutine Tolerance_For_Gauss_Gradients(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real,   intent(out) :: val
  logical,   optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('TOLERANCE_FOR_GAUSS_GRADIENTS',  &
                                 1.0e-3, val, verbose)

  end subroutine
