!==============================================================================!
  recursive subroutine Short_Carry_Short(Sort, a, b)
!------------------------------------------------------------------------------!
!   Quick sort one real array and carry an integer arral along                 !
!                                                                              !
!   Adapted from: https://gist.github.com/1AdAstra1  (good work Olga)          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Sort_Type) :: Sort
  integer(SP)      :: a(:)
  integer(SP)      :: b(:)
!-----------------------------------[Locals]-----------------------------------!
  integer(SP) :: x
  integer     :: i, j, n
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
    call Swap_Short(a(i), a(j))
    call Swap_Short(b(i), b(j))

    i = i + 1
    j = j - 1
  end do

  if(1 < i - 1) call Sort % Short_Carry_Short(a(1:i-1), b(1:i-1))
  if(j + 1 < n) call Sort % Short_Carry_Short(a(j+1:n), b(j+1:n))

  end subroutine
