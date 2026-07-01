!==============================================================================!
  character(SL) function Single_Word_From_Keyboard(File, key_log_entry)
!------------------------------------------------------------------------------!
!>  Reads a single word form keybaord, discarding comments.
!>  A comment is each line which begins with "!", "#" or "%".
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type)           :: File           !! parent class
  character(len=*), optional :: key_log_entry  !! optional for key.log
!==============================================================================!

  ! Remember: 5 is reserved unit for keyboard input
  if(present(key_log_entry)) then
    call File % Read_Line(5, key_log_entry=key_log_entry)
  else
    call File % Read_Line(5)
  end if

  ! That single word is the result of this function
  Single_Word_From_Keyboard = Line % tokens(1)

  end function
