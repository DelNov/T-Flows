!==============================================================================!
  subroutine Normalization_For_Energy_Solver(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real                :: val
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Real_Item('NORMALIZATION_FOR_ENERGY_SOLVER',  &
                                 1.0e-6, val, verbose)

  end subroutine
