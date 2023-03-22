!==============================================================================!
  subroutine Framed(Msg, width, header_text, message_text, one_proc)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Message_Type)           :: Msg
  integer,           intent(in) :: width
  character(*),      intent(in) :: header_text
  character(*),      intent(in) :: message_text
  logical, optional, intent(in) :: one_proc
!-----------------------------------[Locals]-----------------------------------!
  integer       :: w
  character(DL) :: line
!==============================================================================!

  ! Check if all processors have to print
  if(present(one_proc)) then
    if(one_proc) then
      if(this_proc > 1) return
    end if
  end if

  ! Adjust width, if necessary
  w = max(width, len_trim(header_text)+3)

  !------------------------------!
  !   Write the top line first   !
  !------------------------------!
  call Msg % Thick_Line(w)

  !----------------------!
  !   Write the header   !  (you should check it is not too long)
  !----------------------!
  if(header_text .ne. '') then
    line      = ' '
    line(2:2) = '#'
    write(line(4:len_trim(header_text)+5), '(a)')  trim(header_text)
    print '(a)', trim(line)

    call Msg % Dashed_Line(w)
  end if

  !-----------------------------------------------------------!
  !   Write the message text wrapping it into desired width   !
  !-----------------------------------------------------------!
  call Msg % Frameless(w, message_text)

  !--------------------------------!
  !   Write the bottom line last   !
  !--------------------------------!
  call Msg % Thin_Line(w)

  end subroutine
