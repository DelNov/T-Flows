!==============================================================================!
  subroutine Simple_Underrelaxation_For_Momentum(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real,   intent(out) :: val
  logical,   optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('SIMPLE_UNDERRELAXATION_FOR_MOMENTUM', 0.6,  &
                                 val, verbose)

  end subroutine
