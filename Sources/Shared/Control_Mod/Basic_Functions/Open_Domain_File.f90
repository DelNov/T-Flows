!==============================================================================!
  subroutine Open_Domain_File(Control, dom, file_name)
!------------------------------------------------------------------------------!
!   Opens control file for a domain.                                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)             :: Control
  integer,             intent(in) :: dom
  character(len=*),    intent(in) :: file_name
!==============================================================================!

  call File % Open_For_Reading_Ascii(file_name,                     &
                                     Control % dom_file_unit(dom),  &
                                     processor = This_Proc())

  end subroutine
