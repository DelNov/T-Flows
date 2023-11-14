!==============================================================================!
  subroutine Thermal_Conductivity(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads the value of thermal conductivity from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: val      !! value of thermal conductivity
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Real_Item('THERMAL_CONDUCTIVITY', 1.0, val, verbose)

  end subroutine
