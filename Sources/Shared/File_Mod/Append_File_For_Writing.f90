!==============================================================================!
  subroutine File_Mod_Append_File_For_Writing(name_o, file_unit, processor)
!------------------------------------------------------------------------------!
!   Opens file for writing in the first available unit.                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*)  :: name_o
  integer           :: file_unit
  integer, optional :: processor
!-----------------------------------[Locals]-----------------------------------!
  logical       :: file_exists          ! file exists?
  logical, save :: first_call = .true.  ! first call?
!==============================================================================!

  inquire(file = trim(name_o), exist = file_exists)

  open(newunit = file_unit, file = trim(name_o), position = 'append')

  ! All this construct to create the opening message
  if(.not. present(processor)) then
    if(.not. file_exists) then
      print *, '# Creating the file in append mode: ', trim(name_o)
    end if
    if(file_exists .and. first_call) then
      print *, '# Appending to file: ', trim(name_o)
    end if
  else
    if(processor < 2) then
      if(.not. file_exists) then
        print *, '# Creating the file in append mode: ', trim(name_o)
      end if
      if(file_exists .and. first_call) then
        print *, '# Appending to file: ', trim(name_o)
      end if
    end if
  end if

  first_call = .false.

  end subroutine
