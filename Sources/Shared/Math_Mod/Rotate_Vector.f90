!==============================================================================!
  pure function Rotate_Vector(Math, v, k, theta)
!------------------------------------------------------------------------------!
!>  Rotates vector v around the axis k with angle theta
!>  using the Rodrigues' formula.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Math_Type),   intent(in) :: Math           !! parent class
  real, dimension(3)             :: Rotate_Vector  !! result, rotated vector
  real, dimension(3), intent(in) :: v, k           !! input vector
  real,               intent(in) :: theta          !! rotation angle
!-----------------------------------[Locals]-----------------------------------!
  real :: k_dot_v
  real :: k_vec_v(3)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Math)
!==============================================================================!

  k_dot_v = dot_product(k(1:3), v(1:3))
  k_vec_v = Math % Cross_Product(k(1:3), v(1:3))

  Rotate_Vector(1:3) = v(1:3)       * cos(theta)            &
                     + k_vec_v(1:3) * sin(theta)            &
                     + k_dot_v * k(1:3) * (1.0-cos(theta))

  end function
