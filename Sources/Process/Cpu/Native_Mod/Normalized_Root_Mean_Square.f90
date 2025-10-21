!==============================================================================!
  real function Normalized_Root_Mean_Square(Nat, ni, r, A)
!------------------------------------------------------------------------------!
!>  The Normalized_Root_Mean_Square function calculates the root mean square
!>  (RMS) of a given vector 'r', normalized with matrix A diagonal.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type),        intent(in) :: Nat   !! parent class
  integer,                   intent(in) :: ni    !! number of uknowns
  real,                      intent(in) :: r(:)  !! input vector
  type(Matrix_Type), target, intent(in) :: A     !! system matrix
!-----------------------------------[Locals]-----------------------------------!
  real                          :: rms, x_max, x_min, x_max_min
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

  call Global % Sum_Real(rms)

  Normalized_Root_Mean_Square = sqrt(rms)

  end function
