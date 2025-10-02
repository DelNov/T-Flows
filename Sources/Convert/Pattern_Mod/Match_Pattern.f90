!==============================================================================!
  logical function Match_Pattern(Pat, b)
!------------------------------------------------------------------------------!
!>  This subroutine is responsible for comparing a given byte sequence with
!>  the stored pattern in the Pattern_Type instance. It returns a logical value
!>  indicating whether the provided byte sequence matches the pattern.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Pattern_Type)    :: Pat                !! parent Pattern_Type object
  integer(1), intent(in) :: b(1:Pat % length)  !! byte sequence to compare
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

