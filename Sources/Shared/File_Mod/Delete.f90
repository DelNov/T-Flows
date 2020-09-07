!==============================================================================!
  subroutine File_Mod_Delete(name_d)
!------------------------------------------------------------------------------!
!   Delets a file specified by its name.                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*) :: name_d
!-----------------------------------[Locals]-----------------------------------!
  integer :: file_unit
  logical :: file_exists
!==============================================================================!

  ! First check if the file exists
  inquire(file  = trim(name_d),  &
          exist = file_exists)

  ! File exists
  if(file_exists) then
    open(newunit = file_unit, file = trim(name_d), status = 'old')
    close(file_unit, status='delete')
  end if

  end subroutine
