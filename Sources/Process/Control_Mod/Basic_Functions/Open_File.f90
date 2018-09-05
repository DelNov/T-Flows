!==============================================================================!
  subroutine Control_Mod_Open_File
!------------------------------------------------------------------------------!
!   It rarely gets simpler than this - just opens a file for reading.          !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  open(unit   = CONTROL_FILE_UNIT,  &
       file   = CONTROL_FILE_NAME,  &
       action = 'read')

  end subroutine
