!======================================================================!
  real function Distance_Squared(Math, x_a, y_a, z_a, x_b, y_b, z_b)
!----------------------------------------------------------------------!
!  Calculates squared distance between two points in three-dimensions. !
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Arguments]------------------------------!
  class(Math_Type) :: Math
  real             :: x_a, y_a, z_a, x_b, y_b, z_b
!======================================================================!

  Distance_Squared = (x_a-x_b)*(x_a-x_b) &
                   + (y_a-y_b)*(y_a-y_b) &
                   + (z_a-z_b)*(z_a-z_b)

  end function
