!==============================================================================!
  subroutine Warning(Message, width, message_text, file, line, one_proc)
!------------------------------------------------------------------------------!
!>  Displays a warning message in a structured format. It optionally includes
!>  the file name and line number where the warning was triggered.  It is
!>  designed to work in parallel computing environments, allowing for the
!>  option to print from only one processor.
!------------------------------------------------------------------------------!
!   Functionality:                                                             !
!                                                                              !
!   * Header formation:                                                        !
!     - Depending on the presence of file and line arguments, the subroutine   !
!       forms a header string. It concatenates "WARNING" with the file name    !
!       and line number if provided, or just displays "WARNING!" otherwise.    !
!   * Width adjustment:                                                        !
!     - The width (wd) is adjusted to accommodate the length of the header     !
!       text plus a margin.                                                    !
!   * Warning message printing:                                                !
!     - The subroutine prints the warning message using the Framed method of   !
!       the Message object. If one_proc is provided and set to true, it        !
!       checks if the current processor is the first processor (First_Proc()). !
!       If so, or if one_proc is not provided or set to false, it proceeds     !
!       to call Message % Framed.                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Message_Type)           :: Message       !! parent class
  integer, intent(in)           :: width         !! width of the message
  character(*)                  :: message_text  !! actual message thext
  character(*),        optional :: file          !! file with warning
  integer, intent(in), optional :: line          !! line with warning
  logical,             optional :: one_proc      !! print from one processor?
!-----------------------------------[Locals]-----------------------------------!
  integer       :: wd
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


  end subroutine
