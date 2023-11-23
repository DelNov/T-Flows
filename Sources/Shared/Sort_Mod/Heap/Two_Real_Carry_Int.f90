!==============================================================================!
  pure subroutine Two_Real_Carry_Int(Sort, a1, a2, b)
!------------------------------------------------------------------------------!
!>  Heap sort two real arrays and carry one integer array along.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Sort_Type), intent(inout) :: Sort          !! parent class
  real,             intent(inout) :: a1(:), a2(:)  !! array for sorting
  integer,          intent(inout) :: b(:)          !! array to carry on
!-----------------------------------[Locals]-----------------------------------!
  integer :: n, i, ir, j, l
  real    :: a_1, a_2
  integer :: b_
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
      b_  = b (l)
    else
      a_1 = a1(ir)
      a_2 = a2(ir)
      b_  = b (ir)
      a1(ir) = a1(1)
      a2(ir) = a2(1)
      b (ir) = b (1)
      ir = ir - 1
      if(ir == 1) then
        a1(1) = a_1
        a2(1) = a_2
        b (1) = b_
        return
      end if
    end if

    i = l
    j = 2 * l

    do while(j <= ir)
      if(j < ir) then
        if( Math % Smaller_Real(a1(j), a1(j+1)) .or.   &
            Math % Approx_Real (a1(j), a1(j+1)) .and.  &
            Math % Smaller_Real(a2(j), a2(j+1)) ) then
          j = j + 1
        end if
      end if
      if( Math % Smaller_Real(a_1, a1(j)) .or.   &
          Math % Approx_Real (a_1, a1(j)) .and.  &
          Math % Smaller_Real(a_2, a2(j)) ) then
        a1(i) = a1(j)
        a2(i) = a2(j)
        b (i) = b (j)
        i = j
        j = 2 * j
      else
        exit
      end if
    end do

    a1(i) = a_1
    a2(i) = a_2
    b (i) = b_

  end do

  end subroutine
