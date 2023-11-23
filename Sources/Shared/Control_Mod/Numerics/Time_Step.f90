!==============================================================================!
  subroutine Time_Step(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads the time step from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: val      !! value of the time step
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Real_Item('TIME_STEP', 1.0e-2,      &
                                 val, verbose)

  end subroutine
