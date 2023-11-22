!==============================================================================!
  pure subroutine Reverse_Order_Real(Sort, a)
!------------------------------------------------------------------------------!
!   Put a real array in reverse order                                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Sort_Type), intent(in)    :: Sort
  real,             intent(inout) :: a(:)
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j, n
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Sort)
!==============================================================================!

  n = size(a, 1)

  ! Array too small, nothing to do
  if(n < 2) return

  ! Put in reverse order
  do i = 1, n / 2
    j = n - i + 1
    call Swap_Real(a(i), a(j))
  end do

  end subroutine
