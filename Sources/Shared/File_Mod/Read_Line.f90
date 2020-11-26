!==============================================================================!
  subroutine File_Mod_Read_Line(un, reached_end)
!------------------------------------------------------------------------------!
!   Reads a line from a file unit un and discards if it is comment.            !
!   In addition, it breaks the line in tokens (individual strings).            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: un  ! unit
  logical, optional :: reached_end
!-----------------------------------[Locals]-----------------------------------!
  integer      :: i, l
  character(6) :: format = '(a000)'
!==============================================================================!
!   A comment is each line which begins with "!", "#", "$", "%" or "*".        !
!   Input line must not exceed DL characters in length (defined in Const_Mod)  !
!------------------------------------------------------------------------------!

  ! If present, assumed the end of file has not been reached
  if(present(reached_end)) then
    reached_end = .false.
  end if

  !-----------------------------------!
  !  Read the whole line into whole   !
  !-----------------------------------!
  write(format(3:5), '(i3.3)') DL
1 read(un,format, end=2) line % whole

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
