!==============================================================================!
  subroutine Error(Message, width, message_text, file, line, one_proc)
!------------------------------------------------------------------------------!
!>  Displays an error message in a structured format, optionally including the
!>  file name and line number where the error occurred. It's tailored to work
!>  in parallel computing environments, allowing for the option to print from
!>  only one processor. The subroutine also signifies a critical stop in the
!>  execution of the program upon encountering an error.
!------------------------------------------------------------------------------!
!   Functionality:                                                             !
!                                                                              !
!   * Header formation:                                                        !
!     - Depending on the presence of file and line arguments, the subroutine   !
!       forms a header string. It concatenates "ERROR" with the file name and  !
!       line number if provided, helping to pinpoint the location of error.    !
!   * Width adjustment:                                                        !
!     - The width (wd) is adjusted to accommodate the length of the header     !
!       text plus a margin.                                                    !
!   * Error message printing:                                                  !
!     - The subroutine prints the error message using the Framed method of     !
!       the Message object. If one_proc is provided and set to true, it        !
!       checks if the current processor is the first processor (First_Proc()). !
!       If so, or if one_proc is not provided or set to false, it proceeds     !
!       to call Message % Framed.                                              !
!   * Critical program stop:                                                   !
!     - After displaying the error message, the subroutine calls               !
!       Global % End_Parallel and then a stop statement.                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Message_Type)                :: Message       !! parent class
  integer,                intent(in) :: width         !! width of the message
  character(*),           intent(in) :: message_text  !! actuall message text
  character(*), optional, intent(in) :: file          !! file with error
  integer,      optional, intent(in) :: line          !! line with error
  logical,      optional, intent(in) :: one_proc      !! print by 1 processor?
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
