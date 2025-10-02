!==============================================================================!
  pure recursive subroutine Two_Int_Carry_Two_Int(Sort, a1, a2, b, c)
!------------------------------------------------------------------------------!
!>  Quick sort two integer arrays and carry two integer arrays along.
!>  Adapted from: https://gist.github.com/1AdAstra1
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Sort_Type), intent(in)    :: Sort          !! parent class
  integer,          intent(inout) :: a1(:), a2(:)  !! array for sorting
  integer,          intent(inout) :: b(:), c(:)    !! array carried on
!-----------------------------------[Locals]-----------------------------------!
  integer :: x1, x2
  integer :: i, j, n
!==============================================================================!

  n = size(a1, 1)
  x1 = a1( (1+n) / 2 )
  x2 = a2( (1+n) / 2 )
  i = 1
  j = n

  do
    do while ( (a1(i).lt.x1)                        .or.  &
               (a1(i).eq.x1) .and. (a2(i).lt.x2) )
      i = i + 1
    end do
    do while ( (x1.lt.a1(j))                        .or.  &
               (x1.eq.a1(j)) .and. (x2.lt.a2(j)) )
      j = j - 1
    end do
    if(i >= j) exit

    ! Swap values in a, b and c
    call Swap_Int(a1(i), a1(j))
    call Swap_Int(a2(i), a2(j))
    call Swap_Int(b(i),  b(j))
    call Swap_Int(c(i),  c(j))

    i = i + 1
    j = j - 1
  end do

  if(1 < i - 1) call Sort % Two_Int_Carry_Two_Int(a1(1:i-1),  &
                                                  a2(1:i-1),  &
                                                  b (1:i-1),  &
                                                  c (1:i-1))
  if(j + 1 < n) call Sort % Two_Int_Carry_Two_Int(a1(j+1:n),  &
                                                  a2(j+1:n),  &
                                                  b (j+1:n),  &
                                                c (j+1:n))

  end subroutine
