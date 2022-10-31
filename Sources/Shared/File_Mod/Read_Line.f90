!==============================================================================!
  subroutine Read_Line(File, un, reached_end, remove)
!------------------------------------------------------------------------------!
!   Reads a Line from a file unit un and discards if it is comment.            !
!   In addition, it breaks the line in tokens (individual strings).            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type)       :: File
  integer                :: un  ! unit
  logical,      optional :: reached_end
  character(*), optional :: remove
!-----------------------------------[Locals]-----------------------------------!
  integer       :: i, j, n
  integer(1)    :: byte
  character(7)  :: format = '(a0000)'
  character(SL) :: fmtd
!==============================================================================!
!   A comment is each line which begins with "!", "#" or "%".                  !
!   Input line must not exceed MAX_TOKENS*2 characters in length               !
!------------------------------------------------------------------------------!

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
    write(format(3:6), '(i4.4)') MAX_TOKENS*2
    read(un, format, end=2) Line % whole
  else
    Line % whole = ''
    do i = 1, MAX_TOKENS*2
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
