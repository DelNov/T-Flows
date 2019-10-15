!======================================================================!
  real function Math_Mod_Distance(x_a, y_a, z_a, x_b, y_b, z_b)
!----------------------------------------------------------------------!
!  Calculates distance between two points in three-dimensional space.  !
!----------------------------------------------------------------------!
  implicit none
!-----------------------------[Arguments]------------------------------!
  real :: x_a, y_a, z_a, x_b, y_b, z_b
!======================================================================!

  Math_Mod_Distance = sqrt(   (x_a-x_b)*(x_a-x_b) &
                            + (y_a-y_b)*(y_a-y_b) &
                            + (z_a-z_b)*(z_a-z_b) )

  end function
