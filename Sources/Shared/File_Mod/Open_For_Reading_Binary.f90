!==============================================================================!
  subroutine Open_For_Reading_Binary(File, name_i, file_unit, verbose)
!------------------------------------------------------------------------------!
!   Opens file for reading in binary format in first available unit.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type)  :: File
  character(len=*)  :: name_i
  integer           :: file_unit
  logical, optional :: verbose    ! used for tests if file is ascii or binary
!-----------------------------------[Locals]-----------------------------------!
  logical       :: file_exists
  logical       :: verb = .true.
  integer       :: f, l
  character(SL) :: name_l   ! l for the last name
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(File)
!==============================================================================!

  if(present(verbose)) verb = verbose

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
  open(newunit = file_unit,      &
       file    = name_i,         &
       form    = 'unformatted',  &
       access  = 'stream',       &
       action  = 'read')

  ! ... and write a message about it
  if(verb) then
    if(First_Proc()) then
      if(Sequential_Run()) then
        print '(a)', ' # Reading the binary file: ' // trim(name_i)
      else
        name_l = name_i
        f = index(name_i, '/')              + 1
        l = index(name_i, '/', back=.true.) - 1
        if(l-f .eq. 0) write(name_l(f:l), '(i1)')   N_Procs()
        if(l-f .eq. 1) write(name_l(f:l), '(i2.2)') N_Procs()
        if(l-f .eq. 2) write(name_l(f:l), '(i3.3)') N_Procs()
        if(l-f .eq. 3) write(name_l(f:l), '(i4.4)') N_Procs()
        if(l-f .eq. 4) write(name_l(f:l), '(i5.5)') N_Procs()
        print '(a)', ' # Reading the binary files: '  //  &
                     trim(name_i) // ' ... ' // trim(name_l)
      end if
    end if
  end if

  end subroutine
