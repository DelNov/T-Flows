!==============================================================================!
  recursive subroutine Sort_Mod_Int_Carry_Real(a, b)
!------------------------------------------------------------------------------!
!   Quick sort one real array and carry an integer arral along                 !
!                                                                              !
!   Adapted from: https://gist.github.com/1AdAstra1  (good work Olga)          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: a(:)
  real    :: b(:)
!-----------------------------------[Locals]-----------------------------------!
  integer :: x
  integer :: i, j, n
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
    call Swap_Int (a(i), a(j))
    call Swap_Real(b(i), b(j))

    i = i + 1
    j = j - 1
  end do

  if (1 < i - 1) call Sort_Mod_Int_Carry_Real(a(1:i-1), b(1:i-1))
  if (j + 1 < n) call Sort_Mod_Int_Carry_Real(a(j+1:n), b(j+1:n))

  end subroutine 
