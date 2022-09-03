!==============================================================================!
  subroutine Tet_Inertia(Math,                       &
                         x1, y1, z1,                 &
                         x2, y2, z2,                 &
                         x3, y3, z3,                 &
                         x4, y4, z4,                 &
                         ix, iy, iz, ixy, ixz, iyz,  &
                         around_node)
!------------------------------------------------------------------------------!
!   Computes the moment of inertia for a tetrahedron defined with nodes 1 - 4. !
!   around the centrod of tetrahedron.  If optional parameter around_node is   !
!   present, it will shif the center of inertia to the specified node.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Math_Type)  :: Math
  real, intent(in)  :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4
  real, intent(out) :: ix, iy, iz, ixy, ixz, iyz
  integer, optional :: around_node
!-----------------------------------[Locals]-----------------------------------!
  real :: xc, yc, zc  ! centroid of the tetrahedron
  real :: xg, yg, zg  ! shift of the moment of inertia
  real :: vol         ! tetrahedron's volume
!==============================================================================!
!   This is, in essence, Fortran implementation of the Tonon's (2004) paper    !
!                                                                              !
!            Tonon:                Here:                                       !
!                                                                              !
!       |  a  -b' -c' |     |  ix  -ixy -ixz |                                 !
!       | -b'  b  -a' |  =  | -ixy  iy  -iyz |                                 !
!       | -c' -a'  c  |     | -ixz -iyz  iz  |                                 !
!                                                                              !
!------------------------------------------------------------------------------!

  !----------------------------------!
  !   Components of inertia tensor   !
  !----------------------------------!

  ! This is a in Tonon's paper
  ix = y1**2 + y2**2 + y3**2 + y4**2  &
     + z1**2 + z2**2 + z3**2 + z4**2  &
     + y1*y2 + y1*y3 + y1*y4          &
     + y2*y3 + y2*y4                  &
     + y3*y4                          &
     + z1*z2 + z1*z3 + z1*z4          &
     + z2*z3 + z2*z4                  &
     + z3*z4

  ! This is b in Tonon's paper
  iy = x1**2 + x2**2 + x3**2 + x4**2  &
     + z1**2 + z2**2 + z3**2 + z4**2  &
     + x1*x2 + x1*x3 + x1*x4          &
     + x2*x3 + x2*x4                  &
     + x3*x4                          &
     + z1*z2 + z1*z3 + z1*z4          &
     + z2*z3 + z2*z4                  &
     + z3*z4

  ! This is c in Tonon's paper
  iz = x1**2 + x2**2 + x3**2 + x4**2  &
     + y1**2 + y2**2 + y3**2 + y4**2  &
     + x1*x2 + x1*x3 + x1*x4          &
     + x2*x3 + x2*x4                  &
     + x3*x4                          &
     + y1*y2 + y1*y3 + y1*y4          &
     + y2*y3 + y2*y4                  &
     + y3*y4

  ! This is b' in Tonon's paper
  ixy = 2*x1*y1 +   x2*y1 +   x3*y1 +   x4*y1  &
      +   x1*y2 + 2*x2*y2 +   x3*y2 +   x4*y2  &
      +   x1*y3 +   x2*y3 + 2*x3*y3 +   x4*y3  &
      +   x1*y4 +   x2*y4 +   x3*y4 + 2*x4*y4

  ! This is c' in Tonon's paper
  ixz = 2*x1*z1 +   x2*z1 +   x3*z1 +   x4*z1  &
      +   x1*z2 + 2*x2*z2 +   x3*z2 +   x4*z2  &
      +   x1*z3 +   x2*z3 + 2*x3*z3 +   x4*z3  &
      +   x1*z4 +   x2*z4 +   x3*z4 + 2*x4*z4

  ! This is a' in Tonon's paper
  iyz = 2*y1*z1 +   y2*z1 +   y3*z1 +   y4*z1  &
      +   y1*z2 + 2*y2*z2 +   y3*z2 +   y4*z2  &
      +   y1*z3 +   y2*z3 + 2*y3*z3 +   y4*z3  &
      +   y1*z4 +   y2*z4 +   y3*z4 + 2*y4*z4

  ! Final expression for tensor components
  vol = abs(Math % Tet_Volume(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4))

  ix  = ix  * vol /  60.0
  iy  = iy  * vol /  60.0
  iz  = iz  * vol /  60.0
  ixy = ixy * vol / 120.0
  ixz = ixz * vol / 120.0
  iyz = iyz * vol / 120.0

  !------------------------------------------------------------------------!
  !   If around_node is present, you want to shift the moment of inertia   !
  !   to the node. (Comes from parallel axis and parallel plane theorem)   !
  !------------------------------------------------------------------------!
  if(present(around_node)) then

    ! Centroid of the tetrahedron
    xc = 0.25 * (x1 + x2 + x3 + x4)
    yc = 0.25 * (y1 + y2 + y3 + y4)
    zc = 0.25 * (z1 + z2 + z3 + z4)

    ! Shift for the node
    if(around_node .eq. 1) then
      xg = xc - x1
      yg = yc - y1
      zg = zc - z1
    end if
    if(around_node .eq. 2) then
      xg = xc - x2
      yg = yc - y2
      zg = zc - z2
    end if
    if(around_node .eq. 3) then
      xg = xc - x3
      yg = yc - y3
      zg = zc - z3
    end if
    if(around_node .eq. 4) then
      xg = xc - x4
      yg = yc - y4
      zg = zc - z4
    end if

    ix  = ix  + vol * (yg**2 + zg**2)
    iy  = iy  + vol * (xg**2 + zg**2)
    iz  = iz  + vol * (xg**2 + yg**2)
    ixy = ixy + vol * xg * yg
    ixz = ixz + vol * xg * zg
    iyz = iyz + vol * yg * zg

  end if

  end subroutine
