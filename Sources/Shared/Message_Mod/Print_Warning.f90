!==============================================================================!
  subroutine Print_Warning(Msg, width, message_text,  &
                           in_file, at_line)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Message_Type)           :: Msg
  integer, intent(in)           :: width
  character(*)                  :: message_text
  character(*),        optional :: in_file
  integer, intent(in), optional :: at_line
!-----------------------------------[Locals]-----------------------------------!
  type(Tokenizer_Type) :: Tok
  integer              :: w
  character(DL)        :: header_text
!==============================================================================!

  !-------------------------------!
  !   Form the header for error   !
  !-------------------------------!
  if(present(in_file) .and. present(at_line)) then
    write(header_text, '(a,i3)')  "WARNING in file: " // in_file //  &
                                  " at line: ", at_line
  else
    write(header_text, '(a,i3)')  "WARNING!"
  end if

  ! Adjust width, if necessary
  w = max(width, len_trim(header_text)+3)

  !-----------------------------------!
  !   Print the body of the message   !
  !-----------------------------------!
  call Msg % Print_Framed_Text(w, header_text, message_text)

  end subroutine
