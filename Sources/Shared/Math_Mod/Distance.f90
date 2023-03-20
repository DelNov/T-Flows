!==============================================================================!
  pure real function Distance(Math, x_a, y_a, z_a,  &
                                    x_b, y_b, z_b)
!------------------------------------------------------------------------------!
!  Calculates distance between two points in three-dimensional space.          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Math_Type), intent(in) :: Math
  real,             intent(in) :: x_a, y_a, z_a, x_b, y_b, z_b
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Math)
!==============================================================================!

  Distance = sqrt(   (x_a-x_b)*(x_a-x_b) &
                   + (y_a-y_b)*(y_a-y_b) &
                   + (z_a-z_b)*(z_a-z_b) )

  end function
