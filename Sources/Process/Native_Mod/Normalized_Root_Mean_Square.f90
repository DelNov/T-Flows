!==============================================================================!
  real function Normalized_Root_Mean_Square(Native, ni, r, A, x, norm)
!------------------------------------------------------------------------------!
!   Calculates root means square of vector r, normalizing it with entries      !
!   in the system matrix (a), values of unknown (x) and optional norm.         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type) :: Native
  integer            :: ni
  real               :: r(:)  ! this may be only in inner cells
  type(Matrix_Type)  :: A
  real               :: x(:)  ! presumably, this goes to buffer cells
  real, optional     :: norm  ! optional number for normalization
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
  !if( (x_max-x_min) < NANO .and. rms < NANO ) then
  !  rms = PICO
  !else
  !  rms = rms / (x_max - x_min + TINY)
  !end if

  ! rms = rms / (x_max - x_min + TINY)
  ! The line above doesn't work
  ! Example: x_max=1.0, x_min=1.0, TINY=1e-30
  !     (x_max - x_min + TINY) = 0.0

  ! New implementation 2023.03.21 Yohei
  ! avoid roundoff error
  !!x_max_min = x_max - x_min
  !!IF (x_max_min==0.0d0) THEN
  !!  x_max_min = x_max_min+TINY
  !!ENDIF
  !!rms = rms / (x_max_min)
  rms = rms / max(x_max-x_min,TINY)

  Normalized_Root_Mean_Square = rms

  end function
