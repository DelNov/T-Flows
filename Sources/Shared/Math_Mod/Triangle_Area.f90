!==============================================================================!
  real function Triangle_Area(Math,     &
                              xa, ya,   &
                              xb, yb,   &
                              xc, yc)
!------------------------------------------------------------------------------!
!>  Area of a tringle in a plane defined with points a, b, c.
!------------------------------------------------------------------------------!
!                                                                              !
!                c                                                             !
!               / \                                                            !
!              /   \                                                           !
!             /     \                                                          !
!            /       \                                                         !
!           /         \                                                        !
!          a-----------b                                                       !
!                                                                              !
!   If the coordinates sent are x and y pairs, it will calculate the area in   !
!   the remaining z direction.  In the same way, if x and z pairs are sent,    !
!   it will return area in y direction.  Finally, if y and z pairs are sent,   !
!   it will return area in x direction.                                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Math_Type) :: Math    !! parent class
  real, intent(in) :: xa, ya  !! point a
  real, intent(in) :: xb, yb  !! point b
  real, intent(in) :: xc, yc  !! point c
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Math)
!==============================================================================!

  Triangle_Area = abs(  xa * (yb-yc)  &
                      + xb * (yc-ya)  &
                      + xc * (ya-yb)  ) * 0.5

  end function
