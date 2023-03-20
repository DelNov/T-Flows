!==============================================================================!
  subroutine Prec_Form(Nat, ni, A, d, d_inv, prec)
!------------------------------------------------------------------------------!
!   Forms preconditioning matrix "D" from provided matrix "A".                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type),        intent(in)  :: Nat
  integer,                   intent(in)  :: ni
  type(Matrix_Type), target, intent(in)  :: A
  real,                      intent(out) :: d(:)
  real,                      intent(out) :: d_inv(:)
  character(SL),             intent(in)  :: prec  ! preconditioner
!-----------------------------------[Locals]-----------------------------------!
  real                          :: sum
  integer                       :: i, j, k
  integer, contiguous,  pointer :: a_col(:), a_row(:), a_dia(:)
  real,    contiguous,  pointer :: a_val(:)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Nat)
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
    !$omp parallel do private(i) shared (d, a_val, a_dia)
    do i = 1, ni
      d(i) = a_val(a_dia(i))
    end do
    !$omp end parallel do

  !--------------------------------------------! 
  !   2) incomplete cholesky preconditioning   !
  !--------------------------------------------!
  else if(prec .eq. 'icc') then
    do i = 1, ni
      sum = a_val(a_dia(i))          ! take diaginal entry
      do j = a_row(i), a_dia(i) - 1  ! only lower traingular
        k = a_col(j)
        sum = sum - d(k) * a_val(j) * a_val(j)
      end do
      d_inv(i) = sum
      d(i)     = 1.0 / sum
    end do

  !---------------------------!
  !   .) no preconditioning   !
  !---------------------------!
  else
    !$omp parallel do private(i) shared (d)
    do i = 1, ni
      d(i) = 1.0
    end do
    !$omp end parallel do
  end if

  end subroutine
