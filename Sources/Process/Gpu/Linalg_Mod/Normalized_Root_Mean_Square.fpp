!==============================================================================!
  real function Normalized_Root_Mean_Square(Lin, n, b, A, c, norm)
!------------------------------------------------------------------------------!
!>  Front-end for calculation sqrt(vector^2 over matrix diagonal^2 operation).
!------------------------------------------------------------------------------!
!   Note: Using intent clause here, was causing slower runs and crashes        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type)        :: Lin   !! parent class
  integer, intent(in)       :: n     !! size of vectors
  real                      :: b(n)  !! vector
  type(Sparse_Type), target :: A     !! matrix
  real                      :: c(n)  !! vector
  real, optional            :: norm  !! normalization factor
!-----------------------------------[Locals]-----------------------------------!
  integer, pointer :: a_dia(:)
  real,    pointer :: a_val(:)
  real             :: c_min, c_max, c_max_min, rms
  integer          :: i, nz
!==============================================================================!

  ! Take aliases
  nz    =  A % nonzeros
  a_dia => A % dia
  a_val => A % val

  rms = 0.0

  !$tf-acc loop begin
  do i = 1, n
    rms = rms + b(i)**2 / a_val(a_dia(i))**2
  end do
  !$tf-acc loop end

  call Global % Sum_Real(rms)

  rms = sqrt(rms)

  if(.not. present(norm)) then
    call Linalg % Min_Max_Vec(n, c_min, c_max, c)
  else
    c_min = 0.0
    c_max = norm
  end if

  call Global % Min_Real(c_min)
  call Global % Max_Real(c_max)

  ! Avoid roundoff error and divided-by-zero
  ! don't do rms = rms / (c_max - c_min + TINY)
  ! because e.g. 1.0 - 1.0 + 1e-30 = 0.0
  c_max_min = c_max - c_min
  c_max_min = max(c_max_min, TINY)

  Normalized_Root_Mean_Square = rms / c_max_min

  end function

