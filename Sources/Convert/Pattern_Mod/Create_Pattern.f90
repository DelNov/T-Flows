!==============================================================================!
  subroutine Create_Pattern(Pat, string)
!------------------------------------------------------------------------------!
!>  The subroutine Create_Pattern is a member of the Pattern_Type and is
!>  responsible for creating a pattern based on a given input string.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Pattern_Type)          :: Pat     !! parent Pattern_Type object
  character(len=*), intent(in) :: string  !! input string
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!==============================================================================!

  Pat % length = len(string)
  allocate(Pat % pattern(Pat % length))

  do i = 1, Pat % length
    Pat % pattern(i) = ichar(string(i:i))
  end do

  end subroutine
