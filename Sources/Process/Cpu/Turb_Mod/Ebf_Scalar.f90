!==============================================================================!
  pure real function Ebf_Scalar(Turb, c, coef)
!------------------------------------------------------------------------------!
!   Calculates blending function for scalars.                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), intent(in) :: Turb
  integer,          intent(in) :: c
  real,             intent(in) :: coef    ! Prandtl number (Shmidt number)
!==============================================================================!

  Ebf_Scalar = 0.01 * ((coef * Turb % y_plus(c)) ** 4  &
             / ((1.0 + 5.0 * coef**3 * Turb % y_plus(c)) + TINY))

  end function
