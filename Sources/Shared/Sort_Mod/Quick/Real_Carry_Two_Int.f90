!==============================================================================!
  pure recursive subroutine Real_Carry_Two_Int(Sort, a, b, c)
!------------------------------------------------------------------------------!
!>  Quick sort one real array and carry two integer arrays along.
!>  Adapted from: https://gist.github.com/1AdAstra1
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Sort_Type), intent(in)    :: Sort        !! parent class
  real,             intent(inout) :: a(:)        !! array for sorting
  integer,          intent(inout) :: b(:), c(:)  !! array to carry on
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

    ! Swap values in a and b
    call Swap_Real(a(i), a(j))
    call Swap_Int (b(i), b(j))
    call Swap_Int (c(i), c(j))

    i = i + 1
    j = j - 1
  end do

  if(1 < i - 1) call Sort % Real_Carry_Two_Int(a(1:i-1), b(1:i-1), c(1:i-1))
  if(j + 1 < n) call Sort % Real_Carry_Two_Int(a(j+1:n), b(j+1:n), c(j+1:n))

  end subroutine
