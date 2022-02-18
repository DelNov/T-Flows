!==============================================================================!
  subroutine Residual_Vector(Native, ni, r, b, A, x)
!------------------------------------------------------------------------------!
!   Calculates residual vector {r} = {b} - [A]{x}                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type) :: Native
  integer            :: ni
  real               :: r(:)  ! this might be only for inner cells
  real               :: b(:)  ! this might be only for inner cells
  type(Matrix_Type)  :: A
  real               :: x(:)  ! this may incude buffer cells
!-----------------------------------[Locals]-----------------------------------!
  integer  :: i, j, k
!==============================================================================!

  !----------------!
  !   r = b - Ax   !
  !----------------!
  ! Why not callig this: call exchange(x) ???
  do i = 1, ni
    r(i) = b(i)
    do j = A % row(i), A % row(i+1) - 1
      k = A % col(j)
      r(i) = r(i) - A % val(j) * x(k)
    end do
  end do

  end subroutine
