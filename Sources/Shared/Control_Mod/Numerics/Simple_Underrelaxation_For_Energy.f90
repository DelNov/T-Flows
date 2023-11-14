!==============================================================================!
  subroutine Simple_Underrelaxation_For_Energy(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads the SIMPLE under-relaxation coefficient for energy (enthalpy).
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: val
  logical,   optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('SIMPLE_UNDERRELAXATION_FOR_ENERGY', 0.5,  &
                                 val, verbose)

  end subroutine
