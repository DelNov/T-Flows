!==============================================================================!
  subroutine Framed(Message, width, header_text, message_text, one_proc)
!------------------------------------------------------------------------------!
!>  Print a formatted message with a header and body, surrounded by framing
!>  lines. The message is formatted to fit within a given width, and various
!>  options are provided for the appearance of the frame and the inclusion of
!>  the header.
!------------------------------------------------------------------------------!
!   Functionality:                                                             !
!   * Processor check (parallel computing context):                            !
!     - If one_proc is provided and set to true, the subroutine checks if the  !
!       current processor is the first processor (First_Proc()) and returns    !
!       without printing if it's not.                                          !
!   * Width adjustment:                                                        !
!     - The width (wd) is adjusted to be the maximum of the provided width or  !
!       the length of the header text plus a margin.                           !
!   * Top line:                                                                !
!     - Prints a thick line at the top of the message frame                    !
!   * Header text:                                                             !
!     - If header_text is not empty, it prints the header within the frame,    !
!       followed by a dashed line. The header is formatted to fit within the   !
!       adjusted width.                                                        !
!   * Message text:                                                            !
!    - The body of the message (message_text) is printed within the frame,     !
!      formatted to wrap within the specified width. This is done using the    !
!      Frameless method of the Message object.                                 !
!   * Bottom line:                                                             !
!     - Finally, a thin line is printed at the bottom of the message frame     !
!       using Message % Thin_Line(wd).                                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Message_Type)           :: Message       !! parent class
  integer,           intent(in) :: width         !! width of the framed message
  character(*),      intent(in) :: header_text   !! text of the header
  character(*),      intent(in) :: message_text  !! main body of the message
  logical, optional, intent(in) :: one_proc      !! printe by one processors?
!-----------------------------------[Locals]-----------------------------------!
  integer       :: wd
  character(DL) :: line
!==============================================================================!

  ! Check if all processors have to print
  if(present(one_proc)) then
    if(one_proc) then
      if(.not. First_Proc()) return
    end if
  end if

  ! Adjust width, if necessary
  wd = max(width, len_trim(header_text)+3)

  !------------------------------!
  !   Write the top line first   !
  !------------------------------!
  call Message % Thick_Line(wd)

  !----------------------!
  !   Write the header   !  (you should check it is not too long)
  !----------------------!
  if(header_text .ne. '') then
    line      = ' '
    line(2:2) = '#'
    write(line(4:len_trim(header_text)+5), '(a)')  trim(header_text)
    print '(a)', trim(line)

    call Message % Dashed_Line(wd)
  end if

  !-----------------------------------------------------------!
  !   Write the message text wrapping it into desired width   !
  !-----------------------------------------------------------!
  call Message % Frameless(wd, message_text)

  !--------------------------------!
  !   Write the bottom line last   !
  !--------------------------------!
  call Message % Thin_Line(wd)

  end subroutine
