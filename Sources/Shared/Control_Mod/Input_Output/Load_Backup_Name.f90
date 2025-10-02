!==============================================================================!
  subroutine Load_Backup_Name(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads the name of the backup file to load from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  character(SL)       :: val      !! name of the backup file to load
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Char_Item('LOAD_BACKUP_NAME', 'skip', val, verbose)

  end subroutine
