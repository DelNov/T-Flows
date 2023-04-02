!==============================================================================!
  subroutine Open_Root_File(Control, file_name)
!------------------------------------------------------------------------------!
!   It rarely gets simpler than this - just opens a file for reading.          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  character(len=*)    :: file_name
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Control)
!==============================================================================!

  call File % Open_For_Reading_Ascii(file_name, root_control_file_unit,  &
                                     processor=This_Proc())

  ! Make root default to begin with
  control_file_unit = root_control_file_unit

  ! Set default values for domain 1
  ! dom_control_file_name(1) = control_file_name
  dom_control_file_unit(1) = root_control_file_unit

  end subroutine
