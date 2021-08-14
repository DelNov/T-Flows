!==============================================================================!
  subroutine Control_Mod_Backup_Save_Interval(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Int_Item('BACKUP_SAVE_INTERVAL', 120, &
                                  val, verbose)

  end subroutine
