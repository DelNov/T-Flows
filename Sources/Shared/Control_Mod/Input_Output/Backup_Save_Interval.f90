!==============================================================================!
  subroutine Backup_Save_Interval(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads backup save interval from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  integer             :: val      !! backup save interval
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Int_Item('BACKUP_SAVE_INTERVAL', 120, val, verbose)

  end subroutine
