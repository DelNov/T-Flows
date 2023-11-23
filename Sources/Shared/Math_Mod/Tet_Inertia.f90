!==============================================================================!
  pure subroutine Tet_Inertia(Math,                             &
                              x_1, y_1, z_1,                    &
                              x_2, y_2, z_2,                    &
                              x_3, y_3, z_3,                    &
                              x_4, y_4, z_4,                    &
                              i_x, i_y, i_z, i_xy, i_xz, i_yz,  &
                              around_node)
!------------------------------------------------------------------------------!
!>  Computes the moment of inertia for a tetrahedron defined with nodes 1 - 4.
!>  around the centrod of tetrahedron.  If optional parameter around_node is
!>  present, it will shif the center of inertia to the specified node.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Math_Type),  intent(in)  :: Math           !! parent class
  real,              intent(in)  :: x_1, y_1, z_1  !! point 1
  real,              intent(in)  :: x_2, y_2, z_2  !! point 2
  real,              intent(in)  :: x_3, y_3, z_3  !! point 3
  real,              intent(in)  :: x_4, y_4, z_4  !! point 4
  real,              intent(out) :: i_x, i_y, i_z, i_xy, i_xz, i_yz
                                          !! tensor of innertia component
  integer, optional, intent(in)  :: around_node
!-----------------------------------[Locals]-----------------------------------!
  real :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4
  real :: xc, yc, zc  ! shift of the coordinates
  real :: vol         ! tetrahedron's volume
!------------------------------------------------------------------------------!
!   This is, in essence, Fortran implementation of the Tonon's (2004) paper    !
!                                                                              !
!            Tonon:                Here:                                       !
!                                                                              !
!       |  a  -b' -c' |     |  i_x  -i_xy -i_xz |                              !
!       | -b'  b  -a' |  =  | -i_xy  i_y  -i_yz |                              !
!       | -c' -a'  c  |     | -i_xz -i_yz  i_z  |                              !
!                                                                              !
!==============================================================================!

  !------------------------------------------------!
  !   Fetch the local coppies of input arguments   !
  !------------------------------------------------!
  x1 = x_1;  x2 = x_2;  x3 = x_3;  x4 = x_4
  y1 = y_1;  y2 = y_2;  y3 = y_3;  y4 = y_4
  z1 = z_1;  z2 = z_2;  z3 = z_3;  z4 = z_4

  !-----------------------------!
  !   Transform local coppies   !
  !-----------------------------!

  ! Argument around_node not present, shift relative to tetrahedron's centroid
  if(.not. present(around_node)) then

    ! Centroid of the tetrahedron
    xc = 0.25 * (x1 + x2 + x3 + x4)
    yc = 0.25 * (y1 + y2 + y3 + y4)
    zc = 0.25 * (z1 + z2 + z3 + z4)

  ! Argument around_node is present, shift coordinates relative to that
  else

    if(around_node .eq. 1) then
      xc = x1
      yc = y1
      zc = z1
    end if
    if(around_node .eq. 2) then
      xc = x2
      yc = y2
      zc = z2
    end if
    if(around_node .eq. 3) then
      xc = x3
      yc = y3
      zc = z3
    end if
    if(around_node .eq. 4) then
      xc = x4
      yc = y4
      zc = z4
    end if

  end if

  x1 = x1 - xc;  x2 = x2 - xc;  x3 = x3 - xc;  x4 = x4 - xc
  y1 = y1 - yc;  y2 = y2 - yc;  y3 = y3 - yc;  y4 = y4 - yc
  z1 = z1 - zc;  z2 = z2 - zc;  z3 = z3 - zc;  z4 = z4 - zc

  !----------------------------------!
  !   Components of inertia tensor   !
  !----------------------------------!

  ! This is a in Tonon's paper
  i_x = y1**2 + y2**2 + y3**2 + y4**2  &
      + z1**2 + z2**2 + z3**2 + z4**2  &
      + y1*y2 + y1*y3 + y1*y4          &
      + y2*y3 + y2*y4                  &
      + y3*y4                          &
      + z1*z2 + z1*z3 + z1*z4          &
      + z2*z3 + z2*z4                  &
      + z3*z4

  ! This is b in Tonon's paper
  i_y = x1**2 + x2**2 + x3**2 + x4**2  &
      + z1**2 + z2**2 + z3**2 + z4**2  &
      + x1*x2 + x1*x3 + x1*x4          &
      + x2*x3 + x2*x4                  &
      + x3*x4                          &
      + z1*z2 + z1*z3 + z1*z4          &
      + z2*z3 + z2*z4                  &
      + z3*z4

  ! This is c in Tonon's paper
  i_z = x1**2 + x2**2 + x3**2 + x4**2  &
      + y1**2 + y2**2 + y3**2 + y4**2  &
      + x1*x2 + x1*x3 + x1*x4          &
      + x2*x3 + x2*x4                  &
      + x3*x4                          &
      + y1*y2 + y1*y3 + y1*y4          &
      + y2*y3 + y2*y4                  &
      + y3*y4

  ! This is b' in Tonon's paper
  i_xy = 2*x1*y1 +   x2*y1 +   x3*y1 +   x4*y1  &
       +   x1*y2 + 2*x2*y2 +   x3*y2 +   x4*y2  &
       +   x1*y3 +   x2*y3 + 2*x3*y3 +   x4*y3  &
       +   x1*y4 +   x2*y4 +   x3*y4 + 2*x4*y4

  ! This is c' in Tonon's paper
  i_xz = 2*x1*z1 +   x2*z1 +   x3*z1 +   x4*z1  &
       +   x1*z2 + 2*x2*z2 +   x3*z2 +   x4*z2  &
       +   x1*z3 +   x2*z3 + 2*x3*z3 +   x4*z3  &
       +   x1*z4 +   x2*z4 +   x3*z4 + 2*x4*z4

  ! This is a' in Tonon's paper
  i_yz = 2*y1*z1 +   y2*z1 +   y3*z1 +   y4*z1  &
       +   y1*z2 + 2*y2*z2 +   y3*z2 +   y4*z2  &
       +   y1*z3 +   y2*z3 + 2*y3*z3 +   y4*z3  &
       +   y1*z4 +   y2*z4 +   y3*z4 + 2*y4*z4

  ! Final expression for tensor components
  vol = abs(Math % Tet_Volume(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4))

  i_x  = i_x  * vol /  60.0 * 6
  i_y  = i_y  * vol /  60.0 * 6
  i_z  = i_z  * vol /  60.0 * 6
  i_xy = i_xy * vol / 120.0 * 6
  i_xz = i_xz * vol / 120.0 * 6
  i_yz = i_yz * vol / 120.0 * 6

  end subroutine
