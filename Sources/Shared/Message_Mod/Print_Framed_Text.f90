!==============================================================================!
  subroutine Print_Framed_Text(Msg, width, header_text, message_text)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Message_Type)      :: Msg
  integer,      intent(in) :: width
  character(*), intent(in) :: header_text
  character(*), intent(in) :: message_text
!-----------------------------------[Locals]-----------------------------------!
  character(DL) :: line
!==============================================================================!

  !------------------------------!
  !   Write the top line first   !
  !------------------------------!
  call Msg % Thick_Line(width)

  !----------------------!
  !   Write the header   !  (you should check it is not too long)
  !----------------------!
  if(header_text .ne. '') then
    line      = ' '
    line(2:2) = '#'
    write(line(4:len_trim(header_text)+5), '(a)')  trim(header_text)
    print '(a)', trim(line)

    call Msg % Dashed_Line(width)
  end if

  !-----------------------------------------------------------!
  !   Write the message text wrapping it into desired width   !
  !-----------------------------------------------------------!
  call Msg % Print_Plain_Text(width, message_text)

  !--------------------------------!
  !   Write the bottom line last   !
  !--------------------------------!
  call Msg % Thin_Line(width)

  end subroutine
