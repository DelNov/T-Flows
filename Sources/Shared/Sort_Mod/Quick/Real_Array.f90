!==============================================================================!
  pure recursive subroutine Real_Array(Sort, a)
!------------------------------------------------------------------------------!
!>  Quick sort one real array.
!>  Adapted from: https://gist.github.com/1AdAstra1
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Sort_Type), intent(in)    :: Sort  !! parent class
  real,             intent(inout) :: a(:)  !! array for sorting
!-----------------------------------[Locals]-----------------------------------!
  real    :: x
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
    call Swap_Real(a(i), a(j))

    i = i + 1
    j = j - 1
  end do

  if(1 < i - 1) call Sort % Real_Array(a(1:i-1))
  if(j + 1 < n) call Sort % Real_Array(a(j+1:n))

  end subroutine
