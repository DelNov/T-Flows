!==============================================================================!
  recursive subroutine Sort_Mod_Short_Carry_Short(a, b)
!------------------------------------------------------------------------------!
!   Quick sort one real array and carry an integer arral along                 !
!                                                                              !
!   Adapted from: https://gist.github.com/1AdAstra1  (good work Olga)          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer(kind=4) :: a(:)
  integer(kind=4) :: b(:)
!-----------------------------------[Locals]-----------------------------------!
  integer(kind=4) :: x
  integer         :: i, j, n
!------------------------------------------------------------------------------!

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
    if (i >= j) exit

    ! Swap values in a and b
    call Swap_Short(a(i), a(j))
    call Swap_Short(b(i), b(j))

    i = i + 1
    j = j - 1
  end do

  if (1 < i - 1) call Sort_Mod_Short_Carry_Short(a(1:i-1), b(1:i-1))
  if (j + 1 < n) call Sort_Mod_Short_Carry_Short(a(j+1:n), b(j+1:n))

  end subroutine
