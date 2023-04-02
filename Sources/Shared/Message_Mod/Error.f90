!==============================================================================!
  subroutine Error(Message, width, message_text, file, line, one_proc)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Message_Type)                :: Message
  integer,                intent(in) :: width
  character(*),           intent(in) :: message_text
  character(*), optional, intent(in) :: file
  integer,      optional, intent(in) :: line
  logical,      optional, intent(in) :: one_proc  ! print from one processor
!-----------------------------------[Locals]-----------------------------------!
  integer       :: wd
  character(DL) :: header_text
!==============================================================================!

  !-------------------------------!
  !   Form the header for error   !
  !-------------------------------!
  if(present(file) .and. present(line)) then
    write(header_text, '(a,i3)')  "ERROR in file: " // file //  &
                                  " at line: ", line
  else if(present(file)) then
    write(header_text, '(a,i3)')  "ERROR in file: " // file
  else
    write(header_text, '(a,i3)')  "ERROR!"
  end if

  ! Adjust width, if necessary
  wd = max(width, len_trim(header_text)+3)

  !-----------------------------------!
  !   Print the body of the message   !
  !-----------------------------------!
  if(present(one_proc)) then
    if(one_proc) then
      if(First_Proc()) call Message % Framed(wd, header_text, message_text)
    else
      call Message % Framed(wd, header_text, message_text)
    end if
  else
    call Message % Framed(wd, header_text, message_text)
  end if

  !----------------------------------------!
  !   Errors are critical by definitiion   !
  !----------------------------------------!
  call Global % End_Parallel
  stop

  end subroutine
