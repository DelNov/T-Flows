!==============================================================================!
  subroutine Open_For_Reading_Binary(File, name_i, file_unit,  &
                                     processor, verbose)
!------------------------------------------------------------------------------!
!   Opens file for reading in binary format in first available unit.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type)  :: File
  character(len=*)  :: name_i
  integer           :: file_unit
  integer, optional :: processor
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  logical :: file_exists
  logical :: verb = .true.
!==============================================================================!

  if(present(verbose)) verb = verbose

  ! First check if the file exists
  inquire(file  = trim(name_i),  &
          exist = file_exists)

  ! File doesn't exist
  if(.not. file_exists) then
    if(.not. present(processor)) then
      call Message % Error(60, "File: " // trim(name_i)    //   &
                               " doesn't exist!"           //   &
                               " This error is critical."  //   &
                               " Exiting!")
    else
      if(processor < 2) then
        call Message % Error(60, "File: " // trim(name_i)    //   &
                                 " doesn't exist!"           //   &
                                 " This error is critical."  //   &
                                 " Exiting!")
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
    if(verb) print *, '# Reading the binary file: ', trim(name_i)
  else
    if(processor < 2) then
      if(verb) print *, '# Reading the binary file: ', trim(name_i)
    end if
  end if

  end subroutine
