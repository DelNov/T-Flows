!==============================================================================!
  subroutine Simple_Underrelaxation_For_Vof(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads the SIMPLE under-relaxation coefficient for VOF.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: val
  logical,   optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('SIMPLE_UNDERRELAXATION_FOR_VOF',  &
                                 0.5, val, verbose)

  end subroutine
