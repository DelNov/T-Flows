!==============================================================================!
  subroutine Read_Line(File, un, reached_end, remove)
!------------------------------------------------------------------------------!
!>  Reads a Line from a file unit un and discards if it is comment and, if
!>  specified so, unwanted characters.  In addition, it breaks the line in
!>  tokens (individual strings).
!>
!>  A comment is each line which begins with "!", "#" or "%".
!>  Input line must not exceed len(Line % whole) characters in length.
!------------------------------------------------------------------------------!
!   Functionality:                                                             !
!                                                                              !
!   * Initial setup:                                                           !
!     - Sets reached_end to .false. if it is provided.                         !
!     - Determines the number of characters to be removed based on remove.     !
!   * Reading the line:                                                        !
!     - Checks if the file is formatted or unformatted using inquire.          !
!     - Reads the line either as a formatted string or character by character, !
!       depending on the file format.                                          !
!     - Ignores carriage return characters (byte .eq. 13) and stops reading    !
!       at newline characters (byte .eq. 10).                                  !
!   * Removing unwanted characters:                                            !
!     - If remove is provided, removes specified characters from the line by   !
!       replacing them with spaces.                                            !
!   * Adjusting the line:                                                      !
!     - Trims leading spaces from the line using adjustl.                      !
!   * Skipping empty or comment lines:                                         !
!     - Skips over empty lines and lines starting with comment characters      !
!   * Tokenization:                                                            !
!     - Calls Line % Parse() to tokenize the line                              !
!   * End of file handling:                                                    !
!     - If the end of the file is reached (indicated by the end=2 label),      !
!       sets reached_end to .true. if it is present.                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type)       :: File         !! parent class
  integer                :: un           !! file unit
  logical,      optional :: reached_end  !! flag to set if the end of the file
                                         !! was reached during the reading
  character(*), optional :: remove       !! list of characters to remove
!-----------------------------------[Locals]-----------------------------------!
  integer       :: i, j, n
  integer(1)    :: byte
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
