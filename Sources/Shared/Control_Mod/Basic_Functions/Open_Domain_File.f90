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

  call File_Mod_Open_File_For_Reading(file_name, dom_control_file_unit(dom))

  end subroutine
