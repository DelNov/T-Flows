!==============================================================================!
  subroutine Residual_Vector(Native, ni, r, b, A, x)
!------------------------------------------------------------------------------!
!   Calculates residual vector {r} = {b} - [A]{x}                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type),         intent(in)  :: Native
  integer,                    intent(in)  :: ni
  real,                       intent(out) :: r(:)  ! only for inner cells
  real,                       intent(in)  :: b(:)  ! only for inner cells
  type(Matrix_Type),  target, intent(in)  :: A
  real,                       intent(in)  :: x(:)  ! may incude buffer cells
!-----------------------------------[Locals]-----------------------------------!
  integer                      :: i, j, k
  real,    contiguous, pointer :: a_val(:)
  integer, contiguous, pointer :: a_col(:), a_row(:)
!==============================================================================!

  ! Take some aliases
  a_col => A % col
  a_row => A % row
  a_val => A % val

  !----------------!
  !   r = b - Ax   !
  !----------------!
  !$omp parallel do private(i, j) shared (r, b, x)
  do i = 1, ni
    r(i) = b(i)
    do j = a_row(i), a_row(i+1) - 1
      k = a_col(j)
      r(i) = r(i) - a_val(j) * x(k)
    end do
  end do
  !$omp end parallel do

  end subroutine
