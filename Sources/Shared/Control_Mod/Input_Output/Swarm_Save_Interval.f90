!==============================================================================!
  subroutine Swarm_Save_Interval(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads results save interval for swarm (particles) from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  integer             :: val      !! swarm (particles) save interval
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Int_Item('SWARM_SAVE_INTERVAL', 60, val, verbose)

  end subroutine
