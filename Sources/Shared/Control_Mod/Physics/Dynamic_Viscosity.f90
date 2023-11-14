!==============================================================================!
  subroutine Dynamic_Viscosity(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads the value of dynamic viscosity from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: val      !! dynamic viscosity value
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Real_Item('DYNAMIC_VISCOSITY', 0.01, val, verbose)

  end subroutine
