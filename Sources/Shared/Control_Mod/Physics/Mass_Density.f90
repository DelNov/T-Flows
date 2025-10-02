!==============================================================================!
  subroutine Mass_Density(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads the value of mass density from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: val      !! mass density value
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Real_Item('MASS_DENSITY', 1.0, val, verbose)

  end subroutine
