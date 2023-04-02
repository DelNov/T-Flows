!==============================================================================!
  subroutine Normalization_For_Scalars_Solver(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real                :: val
  logical,   optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('NORMALIZATION_FOR_SCALARS_SOLVER',  &
                                 1.0e-6, val, verbose)

  end subroutine
