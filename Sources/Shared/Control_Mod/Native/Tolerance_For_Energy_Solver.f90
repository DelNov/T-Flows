!==============================================================================!
  subroutine Tolerance_For_Energy_Solver(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads linear solver tolerance for energy (enthalpy) equation from the
!>  control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: val      !! tolerance
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Real_Item('TOLERANCE_FOR_ENERGY_SOLVER',  &
                                 1.0e-6, val, verbose)

  end subroutine
