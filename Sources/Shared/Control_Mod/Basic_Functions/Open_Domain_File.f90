!==============================================================================!
  subroutine Control_Mod_Open_Domain_File(dom, file_name)
!------------------------------------------------------------------------------!
!   Opens control file for a domain.                                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer          :: dom
  character(len=*) :: file_name
!==============================================================================!

  call File % Open_For_Reading_Ascii(file_name, dom_control_file_unit(dom),  &
                                     processor=this_proc)

  end subroutine
