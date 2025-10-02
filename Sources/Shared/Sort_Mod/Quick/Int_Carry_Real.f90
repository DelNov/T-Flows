!==============================================================================!
  pure recursive subroutine Int_Carry_Real(Sort, a, b)
!------------------------------------------------------------------------------!
!>  Quick sort one integer array and carry a real array along.
!>  Adapted from: https://gist.github.com/1AdAstra1
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Sort_Type), intent(in)    :: Sort  !! parent class
  integer,          intent(inout) :: a(:)  !! array for sorting
  real,             intent(inout) :: b(:)  !! array to carry on
!-----------------------------------[Locals]-----------------------------------!
  integer :: x
  integer :: i, j, n
!==============================================================================!

  n = size(a, 1)
  x = a( (1+n) / 2 )
  i = 1
  j = n

  do
    do while (a(i) < x)
      i = i + 1
    end do
    do while (x < a(j))
      j = j - 1
    end do
    if(i >= j) exit

    ! Swap values in a and b
    call Swap_Int (a(i), a(j))
    call Swap_Real(b(i), b(j))

    i = i + 1
    j = j - 1
  end do

  if(1 < i - 1) call Sort % Int_Carry_Real(a(1:i-1), b(1:i-1))
  if(j + 1 < n) call Sort % Int_Carry_Real(a(j+1:n), b(j+1:n))

  end subroutine
