!==============================================================================!
  subroutine Turbulent_Prandtl_Number(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads turbulent Prandtl number from the control file.  If not specified
!>  in the control file, it takes the defualt value of 0.9.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: val      !! turbulent Prandtl number
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Real_Item('TURBULENT_PRANDTL_NUMBER', 0.9, val, verbose)

  end subroutine
