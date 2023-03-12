!==============================================================================!
  subroutine Prec_Form(Native, ni, A, D, prec)
!------------------------------------------------------------------------------!
!   Forms preconditioning matrix "D" from provided matrix "A".                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type),        intent(in)    :: Native
  integer,                   intent(in)    :: ni
  type(Matrix_Type), target, intent(in)    :: A
  type(Matrix_Type), target, intent(inout) :: D
  character(SL),             intent(in)    :: prec  ! preconditioner
!-----------------------------------[Locals]-----------------------------------!
  real                          :: sum1
  integer                       :: i, j, k
  real,    contiguous,  pointer :: a_val(:), d_val(:)
  integer, contiguous,  pointer :: a_dia(:), d_dia(:)
!==============================================================================!

  ! Take some aliases
  a_val => A % val
  a_dia => A % dia
  d_val => D % val
  d_dia => D % dia

  !---------------------------------! 
  !   1) diagonal preconditioning   !
  !---------------------------------!
  if(prec .eq. 'jacobi') then
    !$omp parallel do private(i) shared (a_val, a_dia, d_val, d_dia)
    do i = 1, ni
      d_val(d_dia(i)) = a_val(a_dia(i))
    end do
    !$omp end parallel do

  !--------------------------------------------! 
  !   2) incomplete cholesky preconditioning   !
  !--------------------------------------------!
  else if(prec .eq. 'icc') then
    do i = 1, ni
      sum1 = A % val(A % dia(i))       ! take diaginal entry   
      do j = A % row(i), A % dia(i)-1  ! only lower traingular
        k = A % col(j)
        sum1 = sum1 - D % val(D % dia(k)) * A % val(j) * A % val(j)
      end do
      D % val(D % dia(i)) = 1.0 / sum1
    end do

  !---------------------------!
  !   .) no preconditioning   !
  !---------------------------!
  else
    !$omp parallel do private(i) shared (d_val, d_dia)
    do i = 1, ni
      D % val(D % dia(i)) = 1.0
    end do
    !$omp end parallel do
  end if

  end subroutine
