!==============================================================================!
  function Math_Mod_Rotate_Vector(v, k, theta)
!------------------------------------------------------------------------------!
!   Rotates vector "v" around the axis "k" with angle "theta"                  !
!   using the Rodrigues' formula.                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, dimension(3)             :: Math_Mod_Rotate_Vector
  real, dimension(3), intent(in) :: v, k
  real                           :: theta
!-----------------------------------[Locals]-----------------------------------!
  real :: k_dot_v
  real :: k_vec_v(3)
!==============================================================================!

  k_dot_v = dot_product(k(1:3), v(1:3))
  k_vec_v = Math_Mod_Cross_Product(k(1:3), v(1:3))

  Math_Mod_Rotate_Vector(1:3) = v(1:3)       * cos(theta)            &
                              + k_vec_v(1:3) * sin(theta)            &
                              + k_dot_v * k(1:3) * (1.0-cos(theta))

  end function
