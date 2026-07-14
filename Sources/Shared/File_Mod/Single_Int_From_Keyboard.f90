!==============================================================================!
  integer function Single_Int_From_Keyboard(File, key_log_entry)
!------------------------------------------------------------------------------!
!>  Reads a single integer form keybaord, discarding comments.
!>  A comment is each line which begins with "!", "#" or "%".
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type)           :: File           !! parent class
  character(len=*), optional :: key_log_entry  !! optional for key.log
!-----------------------------------[Locals]-----------------------------------!
  integer :: val
!==============================================================================!

  ! Remember: 5 is reserved unit for keyboard input
  if(present(key_log_entry)) then
    call File % Read_Line(5, key_log_entry=key_log_entry)
  else
    call File % Read_Line(5)
  end if
  read(Line % tokens(1), *)  val

  ! That single integer is the result of this function
  Single_Int_From_Keyboard = val

  end function
