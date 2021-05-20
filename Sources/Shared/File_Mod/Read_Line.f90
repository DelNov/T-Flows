!==============================================================================!
  subroutine Read_Line(File, un, reached_end, remove)
!------------------------------------------------------------------------------!
!   Reads a line from a file unit un and discards if it is comment.            !
!   In addition, it breaks the line in tokens (individual strings).            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type)       :: File
  integer                :: un  ! unit
  logical,      optional :: reached_end
  character(*), optional :: remove
!-----------------------------------[Locals]-----------------------------------!
  integer       :: i, j, l, n
  integer(1)    :: byte
  character(6)  :: format = '(a000)'
  character(SL) :: fmtd
!==============================================================================!
!   A comment is each line which begins with "!", "#", "$", "%" or "*".        !
!   Input line must not exceed DL characters in length (defined in Const_Mod)  !
!------------------------------------------------------------------------------!

  ! If present, assumed the end of file has not been reached
  if(present(reached_end)) then
    reached_end = .false.
  end if

  ! Take the number of characters to be removed from the line
  n = 0;  if(present(remove)) n = len_trim(remove)

  !-----------------------------------!
  !  Read the whole line into whole   !
  !-----------------------------------!
  inquire(unit=un, formatted=fmtd)
1 continue
  if(fmtd .eq. 'YES') then
    write(format(3:5), '(i3.3)') DL
    read(un, format, end=2) line % whole
  else
    line % whole = ''
    do i = 1, DL
      read(un, end=2) byte
      if(byte .eq. 10) exit
      if(byte .ne. 13) line % whole(i:i) = char(byte)
    end do
  end if

  ! Remove unwanted characters from it (if n .eq. 0, it won't do it)
  do j = 1, n
    do i = 1, len_trim(line % whole)
      if( line % whole(i:i) .eq. remove(j:j) ) then
        line % whole(i:i) = ' '
      end if
    end do
  end do

  ! Shift the whole line to the left (remove leading spaces)
  line % whole = adjustl(line % whole)

  !--------------------!
  !  Skip empty lines  !
  !--------------------!
  if( line % whole .eq. '' ) goto 1 ! see: man ascii

  !----------------------!
  !  Skip comment lines  !
  !----------------------!
  if( trim(line % whole(1:1)) .eq. '!' .or.               &
      trim(line % whole(1:1)) .eq. '#' .or.               &
      trim(line % whole(1:1)) .eq. '%' ) goto 1

  ! Fetch the first and the last character
  l = len_trim(line % whole)
  line % first = line % whole(1:1)
  line % last  = line % whole(l:l)

  !--------------------------------------!
  !  Parse tokens. This is somehow cool  !
  !--------------------------------------!
  line % n_tokens = 0
  if(line % whole(1:1) >= '!') then
    line % n_tokens = 1
    line % s(1)=1
  end if
  do i=1,DL-2
    if( line % whole(i:  i  ) <  '!' .and.  &
        line % whole(i+1:i+1) >= '!') then
      line % n_tokens = line % n_tokens + 1
      line % s(line % n_tokens) = i+1
    end if
    if( line % whole(i  :i  ) >= '!' .and.  &
        line % whole(i+1:i+1) <  '!') then
      line % e(line % n_tokens) = i
    end if
  end do

  ! Chop them up
  do i = 1, line % n_tokens
    line % tokens(i) = line % whole(line % s(i) : line % e(i))
  end do

  return

  ! Error trap, if here, you reached the end of file
2 if(present(reached_end)) then
    reached_end = .true.
  end if

  end subroutine
