!==============================================================================!
  subroutine Open_For_Reading_Ascii(File, name_i, file_unit)
!------------------------------------------------------------------------------!
!   Opens file for writing in the first available unit.                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type) :: File
  character(len=*) :: name_i
  integer          :: file_unit
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
