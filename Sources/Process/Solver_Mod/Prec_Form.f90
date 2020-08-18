!==============================================================================!
  subroutine Prec_Form(ni, a, d, prec)
!------------------------------------------------------------------------------!
!   Forms preconditioning matrix "d" from provided matrix "a".                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: ni
  type(Matrix_Type) :: a
  type(Matrix_Type) :: d
  character(SL)     :: prec  ! preconditioner
!-----------------------------------[Locals]-----------------------------------!
  real    :: sum1
  integer :: i, j, k
!==============================================================================!

  !---------------------------------! 
  !   1) diagonal preconditioning   !
  !---------------------------------!
  if(prec .eq. 'DIAGONAL') then
    do i = 1, ni
      d % val(d % dia(i)) = a % val(a % dia(i))
    end do

  !--------------------------------------------! 
  !   2) incomplete cholesky preconditioning   !
  !--------------------------------------------!
  else if(prec .eq. 'INCOMPLETE_CHOLESKY') then
    do i = 1, ni
      sum1 = a % val(a % dia(i))       ! take diaginal entry   
      do j = a % row(i), a % dia(i)-1  ! only lower traingular
        k = a % col(j)
        sum1 = sum1 - d % val(d % dia(k)) * a % val(j) * a % val(j)  
      end do
      d % val(d % dia(i)) = 1.0 / sum1
    end do

  !---------------------------!
  !   .) no preconditioning   !
  !---------------------------!
  else
    do i = 1, ni
      d % val(d % dia(i)) = 1.0
    end do
  end if

  end subroutine
