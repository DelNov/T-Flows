!==============================================================================!
  subroutine Warning(Msg, width, message_text, file, line, one_proc)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Message_Type)           :: Msg
  integer, intent(in)           :: width
  character(*)                  :: message_text
  character(*),        optional :: file
  integer, intent(in), optional :: line
  logical,             optional :: one_proc  ! print from one processor only
!-----------------------------------[Locals]-----------------------------------!
  integer       :: w
  character(DL) :: header_text
!==============================================================================!

  !-------------------------------!
  !   Form the header for error   !
  !-------------------------------!
  if(present(file) .and. present(line)) then
    write(header_text, '(a,i3)')  "WARNING in file: " // file //  &
                                  " on line: ", line
  else
    write(header_text, '(a,i3)')  "WARNING!"
  end if

  ! Adjust width, if necessary
  w = max(width, len_trim(header_text)+3)

  !-----------------------------------!
  !   Print the body of the message   !
  !-----------------------------------!
  if(present(one_proc)) then
    if(one_proc) then
      if(this_proc < 2) call Msg % Framed(w, header_text, message_text)
    else
      call Msg % Framed(w, header_text, message_text)
    end if
  else
    call Msg % Framed(w, header_text, message_text)
  end if


  end subroutine
