!==============================================================================!
  subroutine Smagorinsky_Constant(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads Smagorinsky constant from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: val      !! value of the Smagorinsky constant
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Real_Item('SMAGORINSKY_CONSTANT', 0.1, val, verbose)

  end subroutine
