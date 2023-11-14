!==============================================================================!
  subroutine Delete(File, name_d)
!------------------------------------------------------------------------------!
!>  Subroutine to delete a file specified by its name.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type) :: File    !! parent class
  character(len=*) :: name_d  !! name of the file to be deleted
!-----------------------------------[Locals]-----------------------------------!
  integer          :: file_unit
  logical          :: file_exists
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(File)
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
