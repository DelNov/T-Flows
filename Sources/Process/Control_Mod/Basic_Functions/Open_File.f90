!==============================================================================!
  subroutine Control_Mod_Open_File
!------------------------------------------------------------------------------!
!   It rarely gets simpler than this - just opens a file for reading.          !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  call File_Mod_Open_File_For_Reading(CONTROL_FILE_NAME, control_file_unit)

  end subroutine
