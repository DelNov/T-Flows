!==============================================================================!
  subroutine Backup_Save_Interval(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  integer             :: val
  logical,   optional :: verbose
!==============================================================================!

  call Control % Read_Int_Item('BACKUP_SAVE_INTERVAL', 120, val, verbose)

  end subroutine
