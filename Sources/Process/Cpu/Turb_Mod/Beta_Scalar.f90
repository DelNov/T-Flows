!==============================================================================!
  pure real function Beta_Scalar(Turb, coef_l, coef_t)
!------------------------------------------------------------------------------!
!   Calculates elliptic blending function for scalars.                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), intent(in) :: Turb
  real,             intent(in) :: coef_l ! laminar Prandtl (or Schmidt) number
  real,             intent(in) :: coef_t ! turbulent Prandtl (or Schmidt) number
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Turb)
!==============================================================================!

  Beta_Scalar = 9.24 * ((coef_l/coef_t)**0.75 - 1.0)  &
                     * (1.0 + 0.28 * exp(-0.007*coef_l/coef_t))

  end function
