!==============================================================================!
  real function Normalized_Residual(n, nb, mat_a, x, r1) 
!------------------------------------------------------------------------------!
!   Calculates normalized residuals.                                           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Comm_Mod
  use Matrix_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: n, nb
  type(Matrix_Type) :: mat_a
  real              :: x(-nb:n), r1(n)  !  [A]{x}={r1}
!-----------------------------------[Locals]-----------------------------------!
  real    :: error, x_max, x_min
  integer :: i
!==============================================================================!

  ! Compute error normalizing it with main diagonal in the system matrix
  error = 0.0
  do i = 1, n
    error = error + r1(i)**2 / mat_a % val(mat_a % dia(i))**2
  end do  
  call Comm_Mod_Global_Sum_Real(error)
  error = sqrt(error)

  ! Normalize it with absolute values of the unknown
  x_min = minval(x(1:n))
  x_max = maxval(x(1:n))
  call Comm_Mod_Global_Min_Real(x_min)
  call Comm_Mod_Global_Max_Real(x_max)

  ! Create a plateau for very small sources and values
  if( (x_max-x_min) < NANO .and. error < NANO ) then
    error = PICO
  else
    error = error / (x_max - x_min + TINY)
  end if

  Normalized_Residual = error

  end function
