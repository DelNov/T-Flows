!==============================================================================!
  subroutine Open_For_Writing_Binary(File, name_o, file_unit)
!------------------------------------------------------------------------------!
!   Opens file for writing in the first available unit.                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type) :: File
  character(len=*) :: name_o
  integer          :: file_unit
!-----------------------------------[Locals]-----------------------------------!
  integer       :: f, l
  character(SL) :: name_b, name_l   ! b for the base, l for the last name
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(File)
!==============================================================================!

  open(newunit = file_unit,      &
       file    = name_o,         &
       form    = 'unformatted',  &
       status  = 'replace',      &
       access  = 'stream')

  if(First_Proc()) then
    if(Sequential_Run()) then
      print '(a)', ' # Creating the binary file: ' // trim(name_o)
    else
      name_l = name_o
      f = index(name_o, '/')              + 1
      l = index(name_o, '/', back=.true.) - 1
      if(l-f .eq. 0) write(name_l(f:l), '(i1)')   N_Procs()
      if(l-f .eq. 1) write(name_l(f:l), '(i2.2)') N_Procs()
      if(l-f .eq. 2) write(name_l(f:l), '(i3.3)') N_Procs()
      if(l-f .eq. 3) write(name_l(f:l), '(i4.4)') N_Procs()
      if(l-f .eq. 4) write(name_l(f:l), '(i5.5)') N_Procs()
      if(trim(name_o) .eq. trim(name_l)) then  ! for creation of pvtu file
        print '(a)', ' # Creating the binary file: ' // trim(name_o)
      else
        print '(a)', ' # Creating the binary files: '  //  &
                     trim(name_o) // ' ... ' // trim(name_l)
      end if
    end if
  end if

  end subroutine
