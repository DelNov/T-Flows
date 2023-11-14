!==============================================================================!
  subroutine Thick_Line(Msg, w)
!------------------------------------------------------------------------------!
!>  Prints a thick line of a specified width.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Message_Type) :: Msg  !! parent class
  integer, intent(in) :: w    !! line width
!-----------------------------------[Locals]-----------------------------------!
  integer       :: i
  character(DL) :: line
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Msg)
!==============================================================================!

  ! _assert_(w <= DL)

  line = ' '
  do i = 2, w+1
    write(line(i:i), '(a1)')  '='
  end do

  write(line(2:2), '(a1)')  '#'

  print '(a)', line(1:w+1)

  end subroutine

