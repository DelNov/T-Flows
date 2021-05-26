!==============================================================================!
  subroutine Check_Side(Surf, s)
!------------------------------------------------------------------------------!
!  This function calculates radii of inscribed and circumscribed circle        !
!  for a given element (int e)                                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf
  integer                  :: s
!-----------------------------------[Locals]-----------------------------------!
  type(Side_Type), pointer :: side
  type(Elem_Type), pointer :: elem(:)
  integer                  :: a, b, c, d
  integer                  :: i_a, j_a, k_a, i_b, j_b, k_b, sum_a, sum_b
!==============================================================================!

  side => Surf % side(s)
  elem => Surf % elem

  a = side % a
  b = side % b
  c = side % c
  d = side % d

  i_a = elem(side % ea) % v(1)
  j_a = elem(side % ea) % v(2)
  k_a = elem(side % ea) % v(3)

  i_b = elem(side % eb) % v(1)
  j_b = elem(side % eb) % v(2)
  k_b = elem(side % eb) % v(3)

  sum_a = i_a + j_a + k_a - c - d - a
  sum_b = i_b + j_b + k_b - c - d - b

  if(sum_a .ne. 0 .or. sum_b .ne. 0) then
    PRINT *, 'WRONG SIDE!'
    WRITE(100+1, '(A,99I7)') 'sum a, b = "', sum_a, sum_b
    WRITE(100+1, '(A,99I7)') 'abcd :', a, b, c, d
    WRITE(100+1, '(A,99I7)') 'ijk_a:', i_a, j_a, k_a
    WRITE(100+1, '(A,99I7)') 'ijk_b:', i_b, j_b, k_b
    STOP
  end if

  end subroutine
