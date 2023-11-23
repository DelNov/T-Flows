!==============================================================================!
  subroutine Tolerance_For_Gauss_Gradients(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads the tolerance for Gauss-gradient computation from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: val      !! tolerance value
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Real_Item('TOLERANCE_FOR_GAUSS_GRADIENTS',  &
                                 1.0e-3, val, verbose)

  end subroutine
