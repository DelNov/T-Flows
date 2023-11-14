!==============================================================================!
  logical function Is_In_Unix_Format(File, file_name)
!------------------------------------------------------------------------------!
!>  Finds out whether a given file uses Unix-style or Windows-style line
!>  endings.  To identify the line-ending format used in a file. Unix-based
!>  systems typically use a single newline character (\n, ASCII code 10) for
!>  line breaks, whereas Windows systems use a carriage return followed by a
!>  newline character (\r\n, ASCII codes 13 and 10, respectively).  This
!>  functions reads 1024 bytes of a file and counts how many times is char(13)
!>  followed by char(10) and if it is less than five times in 1024 bytes, it
!>  concludes that file is in Unix form.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type) :: File       !! parent class
  character(len=*) :: file_name  !! file name
!-----------------------------------[Locals]-----------------------------------!
  integer    :: file_size, file_unit, cnt_13_10, i
  integer(1) :: byte, next
!==============================================================================!

  inquire(file=file_name, size=file_size)

  ! Open the file in binary mode
  call File % Open_For_Reading_Binary(file_name, file_unit, verbose=.false.)

  ! Read a few lines to guess how are new lines being formed
  ! Sequence 13, 10 is used by Windows for line termination
  cnt_13_10 = 0
  do i = 1, min(1024, file_size)  ! don't read beyond the end of file
    read(file_unit) byte
    if(byte .eq. 13) then
      read(file_unit) next
      if(next .eq. 10) then
        cnt_13_10 = cnt_13_10 + 1
      end if
    end if
  end do

  ! Depending on number of occurences of sequence 13, 10, guess the format
  close(file_unit)
  if(cnt_13_10 .lt. 5) then
    Is_In_Unix_Format = .true.
  else
    Is_In_Unix_Format = .false.
  end if

  end function
