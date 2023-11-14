!==============================================================================!
  character(SL) function Single_Word_From_Keyboard(File)
!------------------------------------------------------------------------------!
!>  Reads a single word form keybaord, discarding comments.
!>  A comment is each line which begins with "!", "#" or "%".
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type) :: File  !! parent class
!==============================================================================!

  call File % Read_Line(5)  ! 5 is reserved unit for keyboard input
  Single_Word_From_Keyboard = Line % tokens(1)

  end function
