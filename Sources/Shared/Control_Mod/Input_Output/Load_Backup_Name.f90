!==============================================================================!
  subroutine Control_Mod_Load_Backup_Name(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL)     :: val
  logical, optional :: verbose
!==============================================================================!

  call Control % Read_Char_Item('LOAD_BACKUP_NAME', 'skip', val, verbose)

  end subroutine
