!==============================================================================!
  subroutine Prec_Solve(Solver, ni, a, d, x, b, prec)
!------------------------------------------------------------------------------!
!   Solves the preconditioning system [d]{x}={b}                               !
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
  class(Solver_Type) :: Solver
  integer            :: ni
  type(Matrix_Type)  :: A
  type(Matrix_Type)  :: D
  real               :: x(:)
  real               :: b(:)
  character(SL)      :: prec  ! preconditioner
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j, k
  real    :: sum1
!==============================================================================!

  !---------------------------------!
  !   1) diagonal preconditioning   !
  !---------------------------------!
  if(prec .eq. 'DIAGONAL') then
    do i = 1, ni
      x(i) = b(i)/d % val(d % dia(i))
    end do

  !--------------------------------------------!
  !   2) incomplete cholesky preconditioning   !
  !--------------------------------------------!
  else if(prec .eq. 'INCOMPLETE_CHOLESKY') then

    ! Forward substitutionn
    do i = 1, ni
      sum1 = b(i)
      do j = A % row(i),A % dia(i)-1     ! only the lower triangular
        k = A % col(j)
        sum1 = sum1 - A % val(j)*x(k)
      end do
      x(i) = sum1 * d % val(d % dia(i))  ! BUG ?
    end do

    do i = 1, ni
      x(i) = x(i) / ( d % val(d % dia(i)) + TINY )
    end do

    ! Backward substitution
    do i = ni, 1, -1
      sum1 = x(i)
      do j = A % dia(i)+1, A % row(i+1)-1        ! upper triangular 
        k = A % col(j)
        if(k <= ni) sum1 = sum1 - A % val(j)*x(k)  ! avoid buffer entries
      end do
      x(i) = sum1* d % val(d % dia(i))           ! BUG ?
    end do

  !---------------------------!
  !   .) no preconditioning   !
  !---------------------------!
  else
    do i = 1, ni
      x(i) = b(i)
    end do
  end if

  end subroutine
