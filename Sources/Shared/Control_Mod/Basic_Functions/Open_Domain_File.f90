!==============================================================================!
  subroutine Open_Domain_File(Control, dom, file_name)
!------------------------------------------------------------------------------!
!   Opens control file for a domain.                                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  integer             :: dom
  character(len=*)    :: file_name
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Control)
!==============================================================================!

  call File % Open_For_Reading_Ascii(file_name,                   &
                                     dom_control_file_unit(dom),  &
                                     processor = This_Proc())

  end subroutine
