!==============================================================================!
  integer function Single_Int_From_Keyboard(File)
!------------------------------------------------------------------------------!
!>  Reads a single integer form keybaord, discarding comments.
!>  A comment is each line which begins with "!", "#" or "%".
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type) :: File  !! parent class
!-----------------------------------[Locals]-----------------------------------!
  integer :: val
!==============================================================================!

  call File % Read_Line(5)  ! 5 is reserved unit for keyboard input
  read(Line % tokens(1), *)  val

  Single_Int_From_Keyboard = val

  end function
