!==============================================================================!
  recursive subroutine Sort_Mod_2_Real_Carry_Int(a1, a2, b)
!------------------------------------------------------------------------------!
!   Quick sort two real arrays and carry one integer array along               !
!                                                                              !
!   Adapted from: https://gist.github.com/1AdAstra1  (good work Olga)          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real    :: a1(:), a2(:)
  integer :: b(:)
!-----------------------------------[Locals]-----------------------------------!
  real    :: x1, x2
  integer :: i, j, n
!------------------------------------------------------------------------------!

  n = size(a1, 1)
  x1 = a1( (1+n) / 2 )
  x2 = a2( (1+n) / 2 )
  i = 1
  j = n

  do
    do while ( (a1(i).lt.x1)                                            .or.  &
               Math_Mod_Approx_Real(a1(i),x1) .and. (a2(i).lt.x2) )
      i = i + 1
    end do
    do while ( (x1.lt.a1(j))                                            .or.  &
               Math_Mod_Approx_Real(x1,a1(j)) .and. (x2.lt.a2(j)) )
      j = j - 1
    end do
    if (i >= j) exit

    ! Swap values in a and b
    call Swap_Real(a1(i), a1(j))
    call Swap_Real(a2(i), a2(j))
    call Swap_Int (b (i), b (j))

    i = i + 1
    j = j - 1
  end do

  if (1 < i - 1) call Sort_Mod_2_Real_Carry_Int(a1(1:i-1),  &
                                                a2(1:i-1),  &
                                                b (1:i-1))
  if (j + 1 < n) call Sort_Mod_2_Real_Carry_Int(a1(j+1:n),  &
                                                a2(j+1:n),  &
                                                b (j+1:n))

  end subroutine 
