!==============================================================================!
  subroutine Control_Mod_Open_File(file_name)
!------------------------------------------------------------------------------!
!   It rarely gets simpler than this - just opens a file for reading.          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*) :: file_name
!==============================================================================!

  control_file_name = file_name

  call File_Mod_Open_File_For_Reading(control_file_name, control_file_unit)

  end subroutine
