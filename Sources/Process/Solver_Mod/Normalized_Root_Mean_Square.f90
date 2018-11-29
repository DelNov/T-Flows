!==============================================================================!
  real function Normalized_Root_Mean_Square(ni, r, a, x, norm)
!------------------------------------------------------------------------------!
!   Calculates root means square of vector r, normalizing it with entries      !
!   in the system matrix (a), values of unknown (x) and optional norm.     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Comm_Mod
  use Matrix_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: ni
  real              :: r(1:)  ! this may be only in inner cells
  type(Matrix_Type) :: a
  real              :: x(1:)  ! presumably, this goes to buffer cells
  real, optional    :: norm   ! optional number for normalization
!-----------------------------------[Locals]-----------------------------------!
  real    :: error, x_max, x_min
  integer :: i
!==============================================================================!

  ! Compute error normalizing it with main diagonal in the system matrix
  error = 0.0
  do i = 1, ni
    error = error + r(i)**2 / a % val(a % dia(i))**2
  end do
  call Comm_Mod_Global_Sum_Real(error)
  error = sqrt(error)

  ! Normalize it with absolute values of the unknown
  if(.not. present(norm)) then
    x_min = minval(x(1:ni))
    x_max = maxval(x(1:ni))
  else
    x_min = 0.0
    x_max = norm
  endif
  call Comm_Mod_Global_Min_Real(x_min)
  call Comm_Mod_Global_Max_Real(x_max)

  ! Create a plateau for very small sources and values
  if( (x_max-x_min) < NANO .and. error < NANO ) then
    error = PICO
  else
    error = error / (x_max - x_min + TINY)
  end if

  Normalized_Root_Mean_Square = error

  end function
