!==============================================================================!
  subroutine Prec_Form(Native, ni, A, D, prec)
!------------------------------------------------------------------------------!
!   Forms preconditioning matrix "D" from provided matrix "A".                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type) :: Native
  integer            :: ni
  type(Matrix_Type)  :: A
  type(Matrix_Type)  :: D
  character(SL)      :: prec  ! preconditioner
!-----------------------------------[Locals]-----------------------------------!
  real    :: sum1
  integer :: i, j, k
!==============================================================================!

  !---------------------------------! 
  !   1) diagonal preconditioning   !
  !---------------------------------!
  if(prec .eq. 'DIAGONAL') then
    do i = 1, ni
      D % val(D % dia(i)) = A % val(A % dia(i))
    end do

  !--------------------------------------------! 
  !   2) incomplete cholesky preconditioning   !
  !--------------------------------------------!
  else if(prec .eq. 'INCOMPLETE_CHOLESKY') then
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
    do i = 1, ni
      D % val(D % dia(i)) = 1.0
    end do
  end if

  end subroutine
