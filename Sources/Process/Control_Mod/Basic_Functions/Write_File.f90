!==============================================================================!
  subroutine Control_Mod_Write_File
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  open(unit   = CONTROL_FILE_UNIT,  &
       file   = CONTROL_FILE_NAME,  &
       action = 'write')

  close(unit = CONTROL_FILE_UNIT)

  end subroutine
