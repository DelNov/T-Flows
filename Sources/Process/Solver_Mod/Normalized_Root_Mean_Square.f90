!==============================================================================!
  real function Normalized_Root_Mean_Square(ni, r, a, x, norm)
!------------------------------------------------------------------------------!
!   Calculates root means square of vector r, normalizing it with entries      !
!   in the system matrix (a), values of unknown (x) and optional norm.         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: ni
  real              :: r(:)  ! this may be only in inner cells
  type(Matrix_Type) :: a
  real              :: x(:)  ! presumably, this goes to buffer cells
  real, optional    :: norm  ! optional number for normalization
!-----------------------------------[Locals]-----------------------------------!
  real    :: rms, x_max, x_min
  integer :: i
!==============================================================================!

  ! Compute rms normalizing it with main diagonal in the system matrix
  rms = 0.0
  do i = 1, ni
    rms = rms + r(i)**2 / A % val(A % dia(i))**2
  end do
  call Comm_Mod_Global_Sum_Real(rms)
  rms = sqrt(rms)

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
  if( (x_max-x_min) < NANO .and. rms < NANO ) then
    rms = PICO
  else
    rms = rms / (x_max - x_min + TINY)
  end if

  Normalized_Root_Mean_Square = rms

  end function
