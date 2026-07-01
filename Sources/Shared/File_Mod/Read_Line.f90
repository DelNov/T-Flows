!==============================================================================!
  subroutine Read_Line(File, un, reached_end, remove, key_log_entry)
!------------------------------------------------------------------------------!
!>  Reads one meaningful line from file unit un into Line % whole.
!>  Empty lines and comment lines are skipped, and optional characters
!>  specified through remove are replaced with blanks.  After a valid
!>  line is found, it is split into tokens and stored in the global Line
!>  object.
!>
!>  Comment lines are lines whose first non-blank character is "!", "#"
!>  or "%".  The input line must not exceed len(Line % whole) characters.
!>
!>  When reading from standard input (unit 5), the consumed line is also
!>  appended to keyboard_input.log.  If key_log_entry is present, it is
!>  written before the consumed line, making the log usable as a commented
!>  replay script.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type)       :: File           !! parent class
  integer                :: un             !! file unit
  logical,      optional :: reached_end    !! flag to set if the end of the file
                                           !! was reached during the reading
  character(*), optional :: remove         !! list of characters to remove
  character(*), optional :: key_log_entry  !! optional headline for key. log
!-----------------------------------[Locals]-----------------------------------!
  integer       :: i, j, n, key_log_unit
  integer(1)    :: byte
  logical       :: first_call = .true.
  character(10) :: time
  character(8)  :: date
  character(7)  :: format = '(a0000)'
  character(SL) :: fmtd
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(File)
!==============================================================================!

  ! If present, assumed the end of file has not been reached
  if(present(reached_end)) then
    reached_end = .false.
  end if

  ! Take the number of characters to be removed from the Line
  n = 0;  if(present(remove)) n = len_trim(remove)

  !------------------------------------!
  !   Read the whole Line into whole   !
  !------------------------------------!
  inquire(unit=un, formatted=fmtd)
1 continue
  if(fmtd .eq. 'YES') then
    write(format(3:6), '(i4.4)') len(Line % whole)
    read(un, format, end=2) Line % whole
  else
    Line % whole = ''
    do i = 1, len(Line % whole)
      read(un, end=2) byte
      if(byte .eq. 10) exit
      if(byte .ne. 13) Line % whole(i:i) = char(byte)
    end do
  end if

  ! Remove unwanted characters from it (if n .eq. 0, it won't do it)
  do j = 1, n
    do i = 1, len_trim(Line % whole)
      if( Line % whole(i:i) .eq. remove(j:j) ) then
        Line % whole(i:i) = ' '
      end if
    end do
  end do

  ! Shift the whole Line to the left (remove leading spaces)
  Line % whole = adjustl(Line % whole)

  !----------------------!
  !   Skip empty lines   !
  !----------------------!
  if( Line % whole .eq. '' ) goto 1 ! see: man ascii

  !------------------------!
  !   Skip comment lines   !
  !------------------------!
  if( trim(Line % whole(1:1)) .eq. '!' .or.               &
      trim(Line % whole(1:1)) .eq. '#' .or.               &
      trim(Line % whole(1:1)) .eq. '%' ) goto 1

  !--------------------------!
  !   Save to keyboard log   !
  !--------------------------!
  if(un .eq. 5) then
    call File % Append_For_Writing_Ascii("keyboard_input.log",  &
                                          key_log_unit,         &
                                          processor = 1)
    if(first_call) then
      call date_and_time(date=date, time=time)

      write(key_log_unit, '(a)')  "#=========================================="
      write(key_log_unit, '(a,a4,"-",a2,"-",a2," ",a2,":",a2,":",a2)')  &
                          "# New log, started on: ",                    &
                          date(1:4), date(5:6), date(7:8),              &
                          time(1:2), time(3:4), time(5:6)
      write(key_log_unit, '(a)')  "#------------------------------------------"
      write(key_log_unit, *)
      first_call = .false.
    end if
    if(present(key_log_entry)) then
      write(key_log_unit, '(a)') key_log_entry
    end if
    write(key_log_unit, *)  trim(Line % whole)
    write(key_log_unit, *)
    close(key_log_unit)
  end if

  !----------------------------------------!
  !   Parse tokens. This is somehow cool   !
  !----------------------------------------!
  call Line % Parse()

  return

  !------------------------------------------------------!
  !   Error trap, if here, you reached the end of file   !
  !------------------------------------------------------!
2 if(present(reached_end)) then
    reached_end = .true.
  end if

  end subroutine
