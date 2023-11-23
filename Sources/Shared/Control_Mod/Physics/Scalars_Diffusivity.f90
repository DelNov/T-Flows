!==============================================================================!
  subroutine Scalars_Diffusivity(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads the value of scalar diffusivity from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: val      !! scalar diffusivity value
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Real_Item('SCALARS_DIFFUSIVITY', 1.0e-6, val, verbose)

  end subroutine
