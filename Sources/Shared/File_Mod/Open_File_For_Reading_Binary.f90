!==============================================================================!
  subroutine File_Mod_Open_File_For_Reading_Binary(name_i, file_unit, processor)
!------------------------------------------------------------------------------!
!   Opens file for reading in binary format in first available unit.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*)  :: name_i
  integer           :: file_unit
  integer, optional :: processor
!-----------------------------------[Locals]-----------------------------------!
  logical :: file_exists
!==============================================================================!

  ! First check if the file exists
  inquire(file  = trim(name_i),  &
          exist = file_exists)

  ! File doesn't exist
  if(.not. file_exists) then
    if(.not. present(processor)) then
      print *, '#============================================='
      print *, '# File: "', trim(name_i), '" doesn''t exist!'
      print *, '# This error is critical. Exiting!'
      print *, '#---------------------------------------------'
      stop
    else
      if(processor < 2) then
        print *, '#============================================='
        print *, '# File: "', trim(name_i), '" doesn''t exist!'
        print *, '# This error is critical. Exiting!'
        print *, '#---------------------------------------------'
        stop
      end if
    end if
  end if

  ! File exists; open it ...
  open(newunit = file_unit,      &
       file    = name_i,         &
       form    = 'unformatted',  &
       access  = 'stream',       &
       action  = 'read')

  ! ... and write a message about it
  if(.not. present(processor)) then
    print *, '# Reading the file: ', trim(name_i)
  else
    if(processor < 2) then
      print *, '# Reading the file: ', trim(name_i)
    end if
  end if

  end subroutine
