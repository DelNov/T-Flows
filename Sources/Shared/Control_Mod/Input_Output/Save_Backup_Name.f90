!==============================================================================!
  subroutine Save_Backup_Name(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads the name of the backup file to save from the control file.
!>  (I believe this is not used at all.)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  character(SL)       :: val      !! backup file to save
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Char_Item('SAVE_BACKUP_NAME', 'skip', val, verbose)

  end subroutine
