!==============================================================================!
  subroutine Framed(Message, width, header_text, message_text, one_proc)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Message_Type)           :: Message
  integer,           intent(in) :: width
  character(*),      intent(in) :: header_text
  character(*),      intent(in) :: message_text
  logical, optional, intent(in) :: one_proc
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
