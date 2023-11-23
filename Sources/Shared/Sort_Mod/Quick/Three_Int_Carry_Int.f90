!==============================================================================!
  pure recursive subroutine Three_Int_Carry_Int(Sort, a1, a2, a3, b)
!------------------------------------------------------------------------------!
!>  Quick sort three integer arrays and carry an integer array along.
!>  Adapted from: https://gist.github.com/1AdAstra1
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Sort_Type), intent(inout) :: Sort                 !! parent class
  integer,          intent(inout) :: a1(:), a2(:), a3(:)  !! array for sorting
  integer,          intent(inout) :: b(:)                 !! array carried on
!-----------------------------------[Locals]-----------------------------------!
  integer :: x1, x2, x3
  integer :: i, j, n
!==============================================================================!

  n = size(a1, 1)
  x1 = a1( (1+n) / 2 )
  x2 = a2( (1+n) / 2 )
  x3 = a3( (1+n) / 2 )
  i = 1
  j = n

  do
    do while ( (a1(i).lt.x1)                                            .or.  &
               (a1(i).eq.x1) .and. (a2(i).lt.x2)                        .or.  &
               (a1(i).eq.x1) .and. (a2(i).eq.x2) .and. (a3(i).lt.x3)  )
      i = i + 1
    end do
    do while ( (x1.lt.a1(j))                                            .or.  &
               (x1.eq.a1(j)) .and. (x2.lt.a2(j))                        .or.  &
               (x1.eq.a1(j)) .and. (x2.eq.a2(j)) .and. (x3.lt.a3(j))  )
      j = j - 1
    end do
    if(i >= j) exit

    ! Swap values in a and b
    call Swap_Int(a1(i), a1(j))
    call Swap_Int(a2(i), a2(j))
    call Swap_Int(a3(i), a3(j))
    call Swap_Int(b(i),  b(j))

    i = i + 1
    j = j - 1
  end do

  if(1 < i - 1) call Sort % Three_Int_Carry_Int(a1(1:i-1),  &
                                                a2(1:i-1),  &
                                                a3(1:i-1),  &
                                                b (1:i-1))
  if(j + 1 < n) call Sort % Three_Int_Carry_Int(a1(j+1:n),  &
                                                a2(j+1:n),  &
                                                a3(j+1:n),  &
                                                b (j+1:n))

  end subroutine
