!==============================================================================!
  subroutine Residual_Vector(ni, r, b, a, x)
!------------------------------------------------------------------------------!
!   Calculates residual vector {r} = {b} - [A]{x}                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: ni
  real              :: r(:)  ! this might be only for inner cells
  real              :: b(:)  ! this might be only for inner cells
  type(Matrix_Type) :: a
  real              :: x(:)  ! this may incude buffer cells
!-----------------------------------[Locals]-----------------------------------!
  integer  :: i, j, k
!==============================================================================!

  !----------------!
  !   r = b - Ax   !
  !----------------!
  ! Why not callig this: call exchange(x) ???
  do i = 1, ni
    r(i) = b(i)
    do j = a % row(i), a % row(i+1) - 1
      k = a % col(j)
      r(i) = r(i) - a % val(j) * x(k)
    end do
  end do

  end subroutine
