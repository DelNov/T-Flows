!==============================================================================!
  pure real function Tet_Volume(Math,        &
                                x1, y1, z1,  &
                                x2, y2, z2,  &
                                x3, y3, z3,  &
                                x4, y4, z4)
!------------------------------------------------------------------------------!
!>  Returns the volume of tethraedra spanned with nodes 1 to 4.
!------------------------------------------------------------------------------!
!   The order of nodes matters here, so you should either be very careful to   !
!   send the nodes in the right order (which is a bit of a nuissance) or       !
!   use the absolute value of the volume computed here (that is mostly done).  !
!                                                                              !
!                4-----3                                                       !
!               / \  . |                                                       !
!              /   \   |                                                       !
!             /  .  \  |                                                       !
!            / .     \ |                                                       !
!           /.        \|                                                       !
!          1-----------2                                                       !
!                                                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Math_Type), intent(in) :: Math        !! parent class
  real,             intent(in) :: x1, y1, z1  !! point 1
  real,             intent(in) :: x2, y2, z2  !! point 2
  real,             intent(in) :: x3, y3, z3  !! point 3
  real,             intent(in) :: x4, y4, z4  !! point 4
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Math)
!==============================================================================!

  Tet_Volume = ( ( (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1) ) * (x4-x1) +    &
                 ( (x3-x1)*(z2-z1) - (x2-x1)*(z3-z1) ) * (y4-y1) +    &
                 ( (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1) ) * (z4-z1) )    &
             / 6.0

  end function
