!==============================================================================!
  real function Math_Mod_Tet_Volume(xa, ya, za,  &
                                    xb, yb, zb,  &
                                    xc, yc, zc,  &
                                    xd, yd, zd)
!------------------------------------------------------------------------------!
!   Returns the volume of tethraedra spanned with nodes "a", "b", "c" and "d". !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(in) :: xa, ya, za, xb, yb, zb, xc, yc, zc, xd, yd, zd
!==============================================================================!
!                                                                              !
!   The order of nodes (a,b,c and d) matters.                                  !
!                                                                              !
!                d-----c                                                       !
!               / \  . |                                                       !
!              /   \   |                                                       !
!             /  .  \  |    I am not 100% sure that the figure is OK           !
!            / .     \ |                                                       !
!           /.        \|                                                       !
!          a-----------b                                                       !
!                                                                              !
!------------------------------------------------------------------------------!

  Math_Mod_Tet_Volume = ( ( (yb-ya)*(zc-za) - (yc-ya)*(zb-za) ) * (xd-xa) +    &
                          ( (xc-xa)*(zb-za) - (xb-xa)*(zc-za) ) * (yd-ya) +    &
                          ( (xb-xa)*(yc-ya) - (xc-xa)*(yb-ya) ) * (zd-za) )    &
                      / 6.0

  end function
