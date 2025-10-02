!==============================================================================!
  subroutine Heat_Capacity(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads the value of thermal conductivity from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: val      !! heat capacity
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Real_Item('HEAT_CAPACITY', 1.0, val, verbose)

  end subroutine
