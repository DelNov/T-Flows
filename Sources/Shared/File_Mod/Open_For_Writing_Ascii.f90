!==============================================================================!
  subroutine Open_For_Writing_Ascii(File, name_o, file_unit)
!------------------------------------------------------------------------------!
!   Opens file for writing in the first available unit.                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type) :: File
  character(len=*) :: name_o
  integer          :: file_unit
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(File)
!==============================================================================!

  open(newunit = file_unit, file = trim(name_o), status = 'replace')

  if(First_Proc()) then
    print '(a)', ' # Creating the ASCII file: ' // trim(name_o)
  end if

  end subroutine
