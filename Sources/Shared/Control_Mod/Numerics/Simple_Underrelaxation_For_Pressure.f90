!==============================================================================!
  subroutine Simple_Underrelaxation_For_Pressure(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads the SIMPLE under-relaxation coefficient for pressure.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: val
  logical,   optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('SIMPLE_UNDERRELAXATION_FOR_PRESSURE', 0.4,  &
                                 val, verbose)

  end subroutine
