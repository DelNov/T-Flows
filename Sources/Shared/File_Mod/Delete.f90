!==============================================================================!
  subroutine File_Mod_Delete(name_d)
!------------------------------------------------------------------------------!
!   Delets a file specified by its name.                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*) :: name_d
  integer          :: file_unit
!==============================================================================!

  open(newunit = file_unit, file = trim(name_d), status = 'old')
  close(file_unit, status='delete')

  end subroutine
