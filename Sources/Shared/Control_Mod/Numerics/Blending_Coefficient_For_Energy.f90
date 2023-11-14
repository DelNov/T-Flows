!==============================================================================!
  subroutine Blending_Coefficient_For_Energy(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads blending coefficient for energy (enthalpy) from control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: val      !! blending coefficient value
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Real_Item('BLENDING_COEFFICIENT_FOR_ENERGY', 1.0,  &
                                 val, verbose)

  end subroutine
