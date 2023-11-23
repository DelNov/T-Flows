!==============================================================================!
  subroutine Tolerance_For_Simple_Algorithm(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads the tolerance for SIMPLE algorithm from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: val      !! tolerance value
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Real_Item('TOLERANCE_FOR_SIMPLE_ALGORITHM',  &
                                 1.0e-4, val, verbose)

  end subroutine
