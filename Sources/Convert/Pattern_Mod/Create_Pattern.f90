!==============================================================================!
  subroutine Create_Pattern(Pat, string)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Pattern_Type)          :: Pat
  character(len=*), intent(in) :: string
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!==============================================================================!

  Pat % length = len(string)
  allocate(Pat % pattern(Pat % length))

  do i = 1, Pat % length
    Pat % pattern(i) = ichar(string(i:i))
  end do

  end subroutine
