!==============================================================================!
  subroutine Surface_Tension(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads the value of surface tension coefficient from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: val      !! surface tension value
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Real_Item('SURFACE_TENSION', 0.0, val, verbose)

  end subroutine
