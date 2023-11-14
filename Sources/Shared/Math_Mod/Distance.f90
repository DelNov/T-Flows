!==============================================================================!
  pure real function Distance(Math, x_a, y_a, z_a,  &
                                    x_b, y_b, z_b)
!------------------------------------------------------------------------------!
!>  Calculates distance between two points (a and b) in three-dimensional space.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Math_Type), intent(in) :: Math  !! parent class
  real,             intent(in) :: x_a   !! x coordinate of point a
  real,             intent(in) :: y_a   !! y coordinate of point a
  real,             intent(in) :: z_a   !! z coordinate of point a
  real,             intent(in) :: x_b   !! x coordinate of point b
  real,             intent(in) :: y_b   !! y coordinate of point b
  real,             intent(in) :: z_b   !! z coordinate of point b
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Math)
!==============================================================================!

  Distance = sqrt(   (x_a-x_b)*(x_a-x_b) &
                   + (y_a-y_b)*(y_a-y_b) &
                   + (z_a-z_b)*(z_a-z_b) )

  end function
