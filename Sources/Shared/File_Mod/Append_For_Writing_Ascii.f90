!==============================================================================!
  subroutine Append_For_Writing_Ascii(File, name_o, file_unit, processor)
!------------------------------------------------------------------------------!
!>  Opens a specified file for writing in append mode. If the file doesn't
!>  exist, it's created. The subroutine also prints a message indicating
!>  whether the file is being created or appended to, based on its existence
!>  and other conditions.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type)  :: File       !! parent class
  character(len=*)  :: name_o     !! name of the output file
  integer           :: file_unit  !! file unit to which the file will open
                                  !! dynamically assigned in the function
  integer, optional :: processor  !! processor number
!-----------------------------------[Locals]-----------------------------------!
  logical       :: file_exists          ! file exists?
  logical, save :: first_call = .true.  ! first call?
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(File)
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
