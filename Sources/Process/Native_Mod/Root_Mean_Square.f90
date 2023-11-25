!==============================================================================!
  real function Root_Mean_Square(Nat, ni, r)
!------------------------------------------------------------------------------!
!>  The Root_Mean_Square function calculates the root mean square (RMS) of a
!>  given vector 'r'.  This implementation computes RMS without normalization.
!>  It takes the sum of squares of the vector elements, computes their global
!>  sum in a parallel computing environment, and then takes the square root of
!>  the resulting sum.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type), intent(in) :: Nat   !! parent class
  integer,            intent(in) :: ni    !! number of uknowns (vector length)
  real,               intent(in) :: r(:)  !! input vector for RMS computation
!-----------------------------------[Locals]-----------------------------------!
  real    :: rms
  integer :: i
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Nat)
!==============================================================================!

  ! Compute rms normalizing it with main diagonal in the system matrix
  rms = 0.0
  !$omp parallel do private(i) shared(r) reduction(+ : rms)
  do i = 1, ni
    rms = rms + r(i)**2
  end do
  !$omp end parallel do
  call Global % Sum_Real(rms)
  rms = sqrt(rms)

  Root_Mean_Square = rms

  end function
