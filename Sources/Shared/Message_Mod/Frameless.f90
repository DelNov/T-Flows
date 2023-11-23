!==============================================================================!
  subroutine Frameless(Msg, width, message_text)
!------------------------------------------------------------------------------!
!>  This subroutine formats and prints a message text within a specified width.
!>  It takes a string message and breaks it into lines, ensuring each line does
!>  not exceed the given width.
!------------------------------------------------------------------------------!
!   Functionality:                                                             !
!                                                                              !
!   * Initialization:                                                          !
!     - The message_text is loaded into a tokenizer and parsed.                !
!     - The out_line is initialized to contain a space followed by a # symbol  !
!       at the second position.                                                !
!   * Processing tokens:                                                       !
!     - The subroutine enters a loop over tokens extracted from message_text.  !
!     - Each token is checked to see if it fits in the current line without    !
!       exceeding the specified width.                                         !
!     - If it fits, the token is appended to the out_line.                     !
!     - If it doesn't fit or is a newline character ('\n'), the current line   !
!       is printed, and out_line is reinitialized for the next line.           !
!   * Line management:                                                         !
!     - The line width is dynamically managed by updating cur_p (current       !
!       position) and nex_p (next position) as tokens are added.               !
!     - When a line reaches or exceeds the specified width, it is printed,     !
!       and the process starts for the next line.                              !
!   * Handling new lines:                                                      !
!     - If a newline character is encountered, it causes an immediate line     !
!       break.                                                                 !
!   * Final line print:                                                        !
!     - After processing all tokens, the last line in formation is printed.    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Message_Type)      :: Msg           !! parent class
  integer,      intent(in) :: width         !! maximum width of the text line
  character(*), intent(in) :: message_text  !! the actual text message
!-----------------------------------[Locals]-----------------------------------!
  integer       :: i
  integer       :: cur_p, nex_p
  character(DL) :: out_line
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Msg)
!==============================================================================!

  !---------------------------------------------------!
  !                                                   !
  !   Load the argument into tokenizer and parse it   !
  !                                                   !
  !---------------------------------------------------!
  Line % whole = message_text
  call Line % Parse()

  ! Initialize the output line
  out_line      = ' '
  out_line(2:2) = '#'
  cur_p         = 4    ! starts at 4 to allow for one leading space in output

  !-----------------------------------------------------!
  !                                                     !
  !   Browse through tokens and print them one by one   !
  !                                                     !
  !-----------------------------------------------------!
  i = 1
  do while (i .le. Line % n_tokens)
1   continue
    nex_p = cur_p + len_trim(Line % tokens(i))

    !-------------------------------------------------!
    !   Keep on filling up this line, it still fits   !
    !   (+3 here is to compensate you started at 4)   !
    !-------------------------------------------------!
    if(nex_p < width + 3 .and. Line % tokens(i) .ne. '\n') then
      write(out_line(cur_p:nex_p), '(a)')  trim(Line % tokens(i)) // ' '

    !---------------------------------------------------------!
    !   Line is too long, print it and reset for new tokens   !
    !---------------------------------------------------------!
    else

      ! Print what you have up to now
      print '(a)', trim(out_line)

      ! (Re)initialize the output line
      out_line      = ' '
      out_line(2:2) = '#'
      cur_p         = 4

      ! If you came here because of new line character, skip it.
      if(Line % tokens(i) .eq. '\n') i = i + 1
      goto 1
    end if

    cur_p = nex_p + 1
    i = i + 1
  end do

  !--------------------------------------------------------!
  !                                                        !
  !   Print the last line which was still in the forming   !
  !                                                        !
  !--------------------------------------------------------!
  print '(a)', trim(out_line)

  end subroutine
