!==============================================================================!
  subroutine Check_Side(Surf, s)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to verify the connectivity and integrity of a
!>  side within the surface mesh. Each side is defined by vertices c and d,
!>  and is flanked by two triangular elements ea and eb. The subroutine
!>  identifies the unique vertices a and b in elements ea and eb, respectively,
!>  that are not shared with the side (vertices c and d). It then checks the
!>  sum of the vertex indices for each element to ensure that the side's
!>  connectivity is correct. If the sums sum_a and sum_b are not zero, it
!>  indicates a discrepancy in the mesh structure, triggering an error
!>  message and halting the program to signal potential issues with the mesh
!>  or its construction logic.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf  !! parent class
  integer,      intent(in) :: s     !! side number
!-----------------------------------[Locals]-----------------------------------!
  type(Side_Type), pointer :: side
  type(Elem_Type), pointer :: Elem(:)
  integer                  :: a, b, c, d
  integer                  :: i_a, j_a, k_a, i_b, j_b, k_b, sum_a, sum_b
!==============================================================================!

  side => Surf % side(s)
  Elem => Surf % Elem

  a = side % a
  b = side % b
  c = side % c
  d = side % d

  i_a = Elem(side % ea) % v(1)
  j_a = Elem(side % ea) % v(2)
  k_a = Elem(side % ea) % v(3)

  i_b = Elem(side % eb) % v(1)
  j_b = Elem(side % eb) % v(2)
  k_b = Elem(side % eb) % v(3)

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
