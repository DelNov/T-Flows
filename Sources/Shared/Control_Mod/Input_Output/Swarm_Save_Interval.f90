!==============================================================================!
  subroutine Swarm_Save_Interval(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  integer             :: val
  logical,   optional :: verbose
!==============================================================================!

  call Control % Read_Int_Item('SWARM_SAVE_INTERVAL', 60, val, verbose)

  end subroutine
