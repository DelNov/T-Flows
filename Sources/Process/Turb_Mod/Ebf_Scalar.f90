!==============================================================================!
  real function Ebf_Scalar(Turb, c, coef)
!------------------------------------------------------------------------------!
!   Calculates elliptic blending function for scalars.                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type) :: Turb
  integer          :: c
  real             :: coef    ! Prandtl number (maybe even Shmidt number)
!==============================================================================!

  Ebf_Scalar = 0.01 * ((coef * Turb % y_plus(c)) ** 4  &
             / ((1.0 + 5.0 * coef**3 * Turb % y_plus(c)) + TINY))

  end function
