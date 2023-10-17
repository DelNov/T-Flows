!==============================================================================!
  pure subroutine Three_Int_Carry_Two_Int(Sort, a1, a2, a3, b, c)
!------------------------------------------------------------------------------!
!   Heap sort three integer arrays and carry two integer arrays along          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Sort_Type), intent(inout) :: Sort
  integer,          intent(inout) :: a1(:), a2(:), a3(:)
  integer,          intent(inout) :: b (:), c (:)
!-----------------------------------[Locals]-----------------------------------!
  integer :: n, i, ir, j, l
  integer :: a_1, a_2, a_3
  integer :: b_, c_
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Sort)
!==============================================================================!

  n = size(a1, 1)

  if(n < 2) return

  l = n / 2 + 1
  ir = n

  do
    if(l > 1) then
      l = l - 1
      a_1 = a1(l)
      a_2 = a2(l)
      a_3 = a3(l)
      b_  = b (l)
      c_  = c (l)
    else
      a_1 = a1(ir)
      a_2 = a2(ir)
      a_3 = a3(ir)
      b_  = b (ir)
      c_  = c (ir)
      a1(ir) = a1(1)
      a2(ir) = a2(1)
      a3(ir) = a3(1)
      b (ir) = b (1)
      c (ir) = c (1)
      ir = ir - 1
      if(ir == 1) then
        a1(1) = a_1
        a2(1) = a_2
        a3(1) = a_3
        b (1) = b_
        c (1) = c_
        return
      end if
    end if

    i = l
    j = 2 * l

    do while(j <= ir)
      if(j < ir) then
        if(                                              a1(j) < a1(j+1) .or.  &
           a1(j) == a1(j+1) .and.                        a2(j) < a2(j+1) .or.  &
           a1(j) == a1(j+1) .and. a2(j) == a2(j+1) .and. a3(j) < a3(j+1) ) then
          j = j + 1
        end if
      end if
      if(                                      a_1 < a1(j) .or.  &
         a_1 == a1(j) .and.                    a_2 < a2(j) .or.  &
         a_1 == a1(j) .and. a_2 == a2(j) .and. a_3 < a3(j)) then
        a1(i) = a1(j)
        a2(i) = a2(j)
        a3(i) = a3(j)
        b (i) = b (j)
        c (i) = c (j)
        i = j
        j = 2 * j
      else
        exit
      end if
    end do

    a1(i) = a_1
    a2(i) = a_2
    a3(i) = a_3
    b (i) = b_
    c (i) = c_

  end do

  end subroutine