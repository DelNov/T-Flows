!==============================================================================!
  subroutine Open_Root_File(Control, file_name)
!------------------------------------------------------------------------------!
!   It rarely gets simpler than this - just opens a file for reading.          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  character(len=*)    :: file_name
!==============================================================================!

  call File % Open_For_Reading_Ascii(file_name, Control % root_file_unit,  &
                                     processor=This_Proc())

  ! Make root default to begin with
  Control % file_unit = Control % root_file_unit

  ! Set default values for domain 1
  Control % dom_file_unit(1) = Control % root_file_unit

  end subroutine
