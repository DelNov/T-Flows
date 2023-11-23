!==============================================================================!
  pure real function Harmonic_Mean(Math, a, b)
!------------------------------------------------------------------------------!
!>  Finds harmonic mean of two input arguments, a and b.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Math_Type), intent(in) :: Math  !! parent class
  real,             intent(in) :: a, b  !! values
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Math)
!==============================================================================!

  Harmonic_Mean = 2.0 / (1.0 / a + 1.0 / b)

  end function
