!==============================================================================!
  subroutine Load_Backup_Name(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  character(SL)       :: val
  logical,   optional :: verbose
!==============================================================================!

  call Control % Read_Char_Item('LOAD_BACKUP_NAME', 'skip', val, verbose)

  end subroutine
