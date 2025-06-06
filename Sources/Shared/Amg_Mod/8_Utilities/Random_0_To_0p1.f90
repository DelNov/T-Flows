!==============================================================================!
  double precision function random_0_to_0p1(amg, s)
!------------------------------------------------------------------------------!
!   Function to create "random" sequence of numbers between 0 and 0.1
!
!   Input:   s      - number between 0 and 0.1
!   Output:  random - number between 0 and 0.1
!            s      - s = random
!
!   This is a clever hack to generate a pseudo-random number in the interval
!   [0,0.1] before modern Fortran introduced built-in random_number().
!   Although it works, it is not statistically random - patterns and
!   correlations are likely and is therefore not suitable for serious
!   stochastic simulations.
!
!   I might replace it with: "call random_number(s); s = 0.1d0 * s" one day
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(amg_type)  :: amg
  double precision :: s
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  random_0_to_0p1 = 100.0d0*dexp(s)
  random_0_to_0p1 = random_0_to_0p1 - dble(int(random_0_to_0p1))
  s = random_0_to_0p1

  end function
