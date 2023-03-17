!==============================================================================!
  subroutine Prec_Solve(Native, ni, A, d, d_inv, x, b, prec)
!------------------------------------------------------------------------------!
!   Solves the preconditioning system [D]{x}={b}                               !
!------------------------------------------------------------------------------!
!   Allows preconditioning of the system by:                                   !
!     1. Diagonal preconditioning                                              !
!     2. Incomplete Cholesky preconditioning                                   !
!                                                                              !
!   The type of precondtioning is chosen by setting the variable prec to 0     !
!   (for no preconditioning), 1 (for diagonal preconditioning) or 2 (for       !
!   incomplete Cholesky preconditioning)                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type),         intent(in)  :: Native
  integer,                    intent(in)  :: ni
  type(Matrix_Type),  target, intent(in)  :: A
  real,                       intent(in)  :: d(:)
  real,                       intent(in)  :: d_inv(:)
  real,                       intent(out) :: x(:)
  real,                       intent(in)  :: b(:)
  character(SL),              intent(in)  :: prec  ! preconditioner
!-----------------------------------[Locals]-----------------------------------!
  integer                       :: i, j, k
  real                          :: sum
  integer, contiguous,  pointer :: a_col(:), a_row(:), a_dia(:)
  real,    contiguous,  pointer :: a_val(:)
!==============================================================================!

  ! Take some aliases
  a_col => A % col
  a_row => A % row
  a_dia => A % dia
  a_val => A % val

  !---------------------------------!
  !   1) diagonal preconditioning   !
  !---------------------------------!
  if(prec .eq. 'jacobi') then
    !$omp parallel do private(i) shared (x, b, d_inv)
    do i = 1, ni
      x(i) = b(i) * d_inv(i)
    end do
    !$omp end parallel do

  !--------------------------------------------!
  !   2) incomplete cholesky preconditioning   !
  !--------------------------------------------!
  else if(prec .eq. 'icc') then

    ! Forward substitutionn
    do i = 1, ni
      sum = b(i)
      do j = a_row(i), a_dia(i) - 1              ! only the lower triangular
        k = a_col(j)
        sum = sum - a_val(j) * x(k)
      end do
      x(i) = sum * d(i)
    end do

    !$omp parallel do private(i) shared (x, d_inv)
    do i = 1, ni
      x(i) = x(i) * d_inv(i)
    end do
    !$omp end parallel do

    ! Backward substitution
    do i = ni, 1, -1
      sum = x(i)
      do j = a_dia(i) + 1, a_row(i+1) - 1        ! upper triangular
        k = a_col(j)
        if(k <= ni) sum = sum - a_val(j) * x(k)  ! avoid buffer entries
      end do
      x(i) = sum * d(i)
    end do

  !---------------------------!
  !   .) no preconditioning   !
  !---------------------------!
  else
    !$omp parallel do private(i) shared (x, b)
    do i = 1, ni
      x(i) = b(i)
    end do
    !$omp end parallel do
  end if

  end subroutine
