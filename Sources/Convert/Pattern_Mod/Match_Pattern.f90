!==============================================================================!
  logical function Match_Pattern(Pat, b)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Pattern_Type)    :: Pat
  integer(1), intent(in) :: b(1:Pat % length)
!-----------------------------------[Locals]-----------------------------------!
  logical :: r
  integer :: i
!==============================================================================!

  r = .true.  ! assume they are the same
  do i = 1, Pat % length
    if(Pat % pattern(i) .ne. b(i)) then
      r = .false.
      exit
    end if
  end do

  Match_Pattern = r

  end function

