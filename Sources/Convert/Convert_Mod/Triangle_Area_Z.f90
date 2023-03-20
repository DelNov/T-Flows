!==============================================================================!
  real function Triangle_Area_Z(Convert,  &
                                xa, ya,   &
                                xb, yb,   &
                                xc, yc)
!------------------------------------------------------------------------------!
!   Area of a tringle in xy plane defined with nodes "a", "b", "c"
!                                                                              !
!                c                                                             !
!               / \                                                            !
!              /   \                                                           !
!             /     \                                                          !
!            /       \                                                         !
!           /         \                                                        !
!          a-----------b                                                       !
!                                                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert
  real, intent(in)    :: xa, ya, xb, yb, xc, yc
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Convert)
!------------------------------------------------------------------------------!

  Triangle_Area_Z = abs(xa * (yb-yc) + xb * (yc-ya) + xc * (ya-yb)) * 0.5

  end function
