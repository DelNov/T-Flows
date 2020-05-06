!==============================================================================!
  subroutine Control_Mod_Open_Root_File(file_name)
!------------------------------------------------------------------------------!
!   It rarely gets simpler than this - just opens a file for reading.          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*) :: file_name
!==============================================================================!

  ! Unused
  ! control_file_name = file_name

  call File_Mod_Open_File_For_Reading(file_name, root_control_file_unit)

  ! Make root default to begin with
  control_file_unit = root_control_file_unit

  ! Set default values for domain 1
  ! dom_control_file_name(1) = control_file_name
  dom_control_file_unit(1) = root_control_file_unit

  end subroutine
