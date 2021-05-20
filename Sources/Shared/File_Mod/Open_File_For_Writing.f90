!==============================================================================!
  subroutine Open_For_Writing_Ascii(File, name_o, file_unit, processor)
!------------------------------------------------------------------------------!
!   Opens file for writing in the first available unit.                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type)  :: File
  character(len=*)  :: name_o
  integer           :: file_unit
  integer, optional :: processor
!==============================================================================!

  open(newunit = file_unit, file = trim(name_o), status = 'replace')

  if(.not. present(processor)) then
    print *, '# Creating the file: ', trim(name_o)
  else
    if(processor < 2) then
      print *, '# Creating the file: ', trim(name_o)
    end if
  end if

  end subroutine
