!==============================================================================!
  subroutine Simple_Underrelaxation_For_Turbulence(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real,   intent(out) :: val
  logical,   optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('SIMPLE_UNDERRELAXATION_FOR_TURBULENCE',  &
                                 0.7, val, verbose)

  end subroutine
