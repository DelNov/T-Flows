!==============================================================================!
  subroutine File_Mod_Open_File_For_Writing_Binary(name_o, file_unit, processor)
!------------------------------------------------------------------------------!
!   Opens file for writing in the first available unit.                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*)  :: name_o
  integer           :: file_unit
  integer, optional :: processor
!==============================================================================!

  open(newunit = file_unit,      &
       file    = name_o,         &
       form    = 'unformatted',  &
       access  = 'stream')

  if(.not. present(processor)) then
    print *, '# Creating the file: ', trim(name_o)
  else
    if(processor < 2) then
      print *, '# Creating the file: ', trim(name_o)
    end if
  end if

  end subroutine
