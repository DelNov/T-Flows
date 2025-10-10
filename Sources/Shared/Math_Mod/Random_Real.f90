!==============================================================================!
  pure subroutine Random_Real(Math, seed, u)
!------------------------------------------------------------------------------!
!>  Generates a real random number in the range 0 to 1, using Park-Miller meth.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Math_Type), intent(in)    :: Math   !! parent class
  integer,          intent(inout) :: seed   !! 32-bit; updated on return
  real,             intent(out)   :: u      !! random real in (0,1)
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: A = 16807       ! Park–Miller constants
  integer, parameter :: M = 2147483647  ! 2^31 - 1
  integer, parameter :: Q = 127773      ! M / A
  integer, parameter :: R = 2836        ! M mod A
!-----------------------------------[Locals]-----------------------------------!
  integer :: k, ns
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Math)
!==============================================================================!

  ! Ensure seed is valid: 1..M-1
  if (seed <= 0 .or. seed >= M) seed = mod(abs(seed), M-1) + 1

  ! Schrage method (avoids 32-bit overflow)
  k  = seed / Q
  ns = A * (seed - k * Q) - R * k
  if (ns <= 0) ns = ns + M

  seed = ns
  u    = real(seed) / real(M)   ! (0,1)

end subroutine

