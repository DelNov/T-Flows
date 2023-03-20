!==============================================================================!
  real function Normalized_Root_Mean_Square(Nat, ni, r, A, x, norm)
!------------------------------------------------------------------------------!
!   Calculates root means square of vector r, normalizing it with entries      !
!   in the system matrix (a), values of unknown (x) and optional norm.         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type),        intent(in) :: Nat
  integer,                   intent(in) :: ni
  real,                      intent(in) :: r(:)  ! this may be in inner cells
  type(Matrix_Type), target, intent(in) :: A
  real,                      intent(in) :: x(:)  ! this goes to buffer cells
  real,            optional, intent(in) :: norm  ! number for normalization
!-----------------------------------[Locals]-----------------------------------!
  real                          :: rms, x_max, x_min
  integer                       :: i
  real,    contiguous,  pointer :: a_val(:)
  integer, contiguous,  pointer :: a_dia(:)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Nat)
!==============================================================================!

  ! Take some aliases
  a_val => A % val
  a_dia => A % dia

  ! Compute rms normalizing it with main diagonal in the system matrix
  rms = 0.0
  !$omp parallel do private(i) shared(r, a_val, a_dia) reduction(+ : rms)
  do i = 1, ni
    rms = rms + r(i)**2 / a_val(a_dia(i))**2
  end do
  !$omp end parallel do

  call Comm_Mod_Global_Sum_Real(rms)
  rms = sqrt(rms)

  ! Normalize it with absolute values of the unknown
  if(.not. present(norm)) then
    x_min = +HUGE
    x_max = -HUGE
    !$omp parallel do private(i) shared(x)  &
    !$omp reduction(max : x_max)            &
    !$omp reduction(min : x_min)
    do i = 1, ni
      x_min = min(x_min, x(i))
      x_max = max(x_max, x(i))
    end do
    !$omp end parallel do
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
