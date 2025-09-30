!==============================================================================!
  real function Random_0_To_0p1(Amg, s)
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
!   I might replace it with: "call random_number(s); s = 0.1 * s" one day
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Arguments]----------------------------------!
  class(Amg_Type) :: Amg
  real            :: s
!------------------------------------[Save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  Random_0_To_0p1 = 100.0 * exp(s)
  Random_0_To_0p1 = Random_0_To_0p1 - real(int(Random_0_To_0p1))
  s = Random_0_To_0p1

  end function
