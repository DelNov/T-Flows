!==============================================================================!
  pure recursive subroutine Int_Array(Sort, a)
!------------------------------------------------------------------------------!
!>  Quick sort one integer array.
!>  Adapted from: https://gist.github.com/1AdAstra1
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Sort_Type), intent(in)    :: Sort  !! parent class
  integer,          intent(inout) :: a(:)  !! array for sorting
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

    ! Swap values in a
    call Swap_Int(a(i), a(j))

    i = i + 1
    j = j - 1
  end do

  if(1 < i - 1) call Sort % Int_Array(a(1:i-1))
  if(j + 1 < n) call Sort % Int_Array(a(j+1:n))

  end subroutine
