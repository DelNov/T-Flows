!==============================================================================!
  subroutine Saturation_Temperature(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads the value of saturation temperature from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: val      !! saturation temperature value
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Real_Item('SATURATION_TEMPERATURE', 100., val, verbose)

  end subroutine
