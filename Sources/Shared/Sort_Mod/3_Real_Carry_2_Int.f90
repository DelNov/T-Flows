!==============================================================================!
  recursive subroutine Sort_Mod_3_Real_Carry_2_Int(a1, a2, a3, b, c)
!------------------------------------------------------------------------------!
!   Quick sort three real arrays (think of three coordinates) and carry two    !
!   integer arrays along                                                       !
!                                                                              !
!   Adapted from: https://gist.github.com/1AdAstra1  (good work Olga)          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real    :: a1(:), a2(:), a3(:)
  integer :: b(:), c(:)
!-----------------------------------[Locals]-----------------------------------!
  real    :: x1, x2, x3
  integer :: i, j, n
!------------------------------------------------------------------------------!

  n = size(a1, 1)
  x1 = a1( (1+n) / 2 )
  x2 = a2( (1+n) / 2 )
  x3 = a3( (1+n) / 2 )
  i = 1
  j = n

  do
    do while ( (a1(i).lt.x1)                                            .or.  &
               Math_Mod_Approx_Real(a1(i),x1) .and. (a2(i).lt.x2)       .or.  &
               Math_Mod_Approx_Real(a1(i),x1) .and.                           &
               Math_Mod_Approx_Real(a2(i),x2) .and. (a3(i).lt.x3)  )
      i = i + 1
    end do
    do while ( (x1.lt.a1(j))                                            .or.  &
               Math_Mod_Approx_Real(x1,a1(j)) .and. (x2.lt.a2(j))       .or.  &
               Math_Mod_Approx_Real(x1,a1(j)) .and.                           &
               Math_Mod_Approx_Real(x2,a2(j)) .and. (x3.lt.a3(j))  )
      j = j - 1
    end do
    if (i >= j) exit

    ! Swap values in a, b and c
    call Swap_Real(a1(i), a1(j))
    call Swap_Real(a2(i), a2(j))
    call Swap_Real(a3(i), a3(j))
    call Swap_Int (b (i), b (j))
    call Swap_Int (c (i), c (j))

    i = i + 1
    j = j - 1
  end do

  if (1 < i - 1) call Sort_Mod_3_Real_Carry_2_Int(a1(1:i-1),  &
                                                  a2(1:i-1),  &
                                                  a3(1:i-1),  &
                                                  b (1:i-1),  &
                                                  c (1:i-1))
  if (j + 1 < n) call Sort_Mod_3_Real_Carry_2_Int(a1(j+1:n),  &
                                                  a2(j+1:n),  &
                                                  a3(j+1:n),  &
                                                  b (j+1:n),  &
                                                  c (j+1:n))

  end subroutine 
