!==============================================================================!
  subroutine Wale_Constant(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads Wale constant from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: val      !! value of the Smagorinsky constant
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Real_Item('WALE_CONSTANT', 0.325, val, verbose)

  end subroutine
