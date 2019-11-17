!==============================================================================!
  subroutine Control_Mod_Write_File
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  open(unit   = control_file_unit,  &
       file   = CONTROL_FILE_NAME,  &
       action = 'write')

  close(unit = control_file_unit)

  end subroutine
