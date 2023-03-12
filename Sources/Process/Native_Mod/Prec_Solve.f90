!==============================================================================!
  subroutine Prec_Solve(Native, ni, A, D, x, b, prec)
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
  type(Matrix_Type),  target, intent(in)  :: D
  real,                       intent(out) :: x(:)
  real,                       intent(in)  :: b(:)
  character(SL),              intent(in)  :: prec  ! preconditioner
!-----------------------------------[Locals]-----------------------------------!
  integer                       :: i, j, k
  real                          :: sum1
  real,    contiguous,  pointer :: d_val(:)
  integer, contiguous,  pointer :: d_dia(:)
!==============================================================================!

  ! Take some aliases
  d_val => D % val
  d_dia => D % dia

  !---------------------------------!
  !   1) diagonal preconditioning   !
  !---------------------------------!
  if(prec .eq. 'jacobi') then
    !$omp parallel do private(i) shared (x, b, d_val, d_dia)
    do i = 1, ni
      x(i) = b(i) / d_val(d_dia(i))
    end do
    !$omp end parallel do

  !--------------------------------------------!
  !   2) incomplete cholesky preconditioning   !
  !--------------------------------------------!
  else if(prec .eq. 'icc') then

    ! Forward substitutionn
    do i = 1, ni
      sum1 = b(i)
      do j = A % row(i),A % dia(i)-1     ! only the lower triangular
        k = A % col(j)
        sum1 = sum1 - A % val(j)*x(k)
      end do
      x(i) = sum1 * D % val(D % dia(i))  ! BUG ?
    end do

    do i = 1, ni
      x(i) = x(i) / ( D % val(D % dia(i)) + TINY )
    end do

    ! Backward substitution
    do i = ni, 1, -1
      sum1 = x(i)
      do j = A % dia(i)+1, A % row(i+1)-1          ! upper triangular
        k = A % col(j)
        if(k <= ni) sum1 = sum1 - A % val(j)*x(k)  ! avoid buffer entries
      end do
      x(i) = sum1* D % val(D % dia(i))             ! BUG ?
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
