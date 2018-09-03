!==============================================================================!
  recursive subroutine Sort_Mod_Real(a)
!------------------------------------------------------------------------------!
!   Quick sort one real array.                                                 !
!                                                                              ! 
!   Adapted from: https://gist.github.com/1AdAstra1  (good work Olga)          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real :: a(:)
!-----------------------------------[Locals]-----------------------------------!
  real    :: x, at
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

    ! Swap values in a
    at   = a(i)
    a(i) = a(j)
    a(j) = at

    i = i + 1
    j = j - 1
  end do
  
  if (1 < i - 1) call Sort_Mod_Real(a(1:i-1))
  if (j + 1 < n) call Sort_Mod_Real(a(j+1:n))

  end subroutine 
