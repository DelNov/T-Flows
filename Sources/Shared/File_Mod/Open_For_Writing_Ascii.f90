!==============================================================================!
  subroutine Open_For_Writing_Ascii(File, name_o, file_unit)
!------------------------------------------------------------------------------!
!>  To open an ASCII file for writing. It also prints a message which file is
!>  being created from one (first) processor.
!>  File unit is assigned dynamically when opening the file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type) :: File       !! parent class
  character(len=*) :: name_o     !! name of the output file
  integer          :: file_unit  !! file unit assigned when opening the file
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(File)
!==============================================================================!

  open(newunit = file_unit, file = trim(name_o), status = 'replace')

  if(First_Proc()) then
    print '(a)', ' # Creating the ASCII file: ' // trim(name_o)
  end if

  end subroutine
