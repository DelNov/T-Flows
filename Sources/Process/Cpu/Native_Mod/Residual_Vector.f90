!==============================================================================!
  subroutine Residual_Vector(Nat, ni, r, b, A, x)
!------------------------------------------------------------------------------!
!>  Residual_Vector is a subroutine for calculating the residual vector for
!>  iterative linear solvers. The residual vector is {r} = {b} - [A]{x},
!>  where {b} is the known vector, [A] is the system matrix, and {x} is the
!>  vector of unknowns. This computation is essential to evaluate the accuracy
!>  of the current iteration.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type),         intent(in)  :: Nat   !! parent class
  integer,                    intent(in)  :: ni    !! number of uknowns
  real,                       intent(out) :: r(:)  !! input vector
  real,                       intent(in)  :: b(:)  !! right hand side vector
  type(Matrix_Type),  target, intent(in)  :: A     !! system matrix
  real,                       intent(in)  :: x(:)  !! solution vector
!-----------------------------------[Locals]-----------------------------------!
  integer                      :: i, j, k
  real,    contiguous, pointer :: a_val(:)
  integer, contiguous, pointer :: a_col(:), a_row(:)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Nat)
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
