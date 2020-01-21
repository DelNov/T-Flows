!==============================================================================!
  real function Turb_Mod_Ebf_Scalar(turb, c, coef)
!------------------------------------------------------------------------------!
!   Calculates elliptic blending function for scalars.                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type) :: turb
  integer         :: c
  real            :: coef    ! Prandtl number (maybe even Shmidt number)
!==============================================================================!

  Turb_Mod_Ebf_Scalar = 0.01 * (coef * turb % y_plus(c) ** 4  &
                        / ((1.0 + 5.0 * coef**3 * turb % y_plus(c)) + TINY))

  end function
