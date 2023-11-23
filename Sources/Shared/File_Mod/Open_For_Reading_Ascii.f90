!==============================================================================!
  subroutine Open_For_Reading_Ascii(File, name_i, file_unit)
!------------------------------------------------------------------------------!
!>  To open an ASCII file for reading. It checks if the file exists before
!>  attempting to open it and reports an error if the file does not exist.
!>  It also prints a message which file is being read from one processor.
!>  File unit is assigned dynamically when opening the file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type) :: File       !! parent class
  character(len=*) :: name_i     !! name of the input file
  integer          :: file_unit  !! file unit assigned when opening the file
!-----------------------------------[Locals]-----------------------------------!
  logical :: file_exists
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(File)
!==============================================================================!

  ! First check if the file exists
  inquire(file  = trim(name_i),  &
          exist = file_exists)

  ! File doesn't exist
  if(.not. file_exists) then
    call Message % Error(60, "File: " // trim(name_i)    //  &
                             " doesn't exist!"           //  &
                             " This error is critical."  //  &
                             " Exiting!",                    &
                             one_proc = .true.)
  end if

  ! File exists; open it ...
  open(newunit = file_unit,     &
       file    = trim(name_i),  &
       action  = 'read')

  ! ... and write a message about it
  if(First_Proc()) then
    print '(a)', ' # Reading the ASCII file: ' // trim(name_i)
  end if

  end subroutine
