!==============================================================================!
  subroutine Surf_Mod_Calculate_Curvatures_From_Edges(surf)
!------------------------------------------------------------------------------!
!   Calculates surface curvatures from edges (four points) and distributes     !
!   the values to elements and nodes.                                          !
!                                                                              !
!   This is the lowest order of curvature calculation in the code.             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Surf_Type), target :: surf
!-----------------------------------[Locals]-----------------------------------!
  real, dimension(4,4)  :: a
  real, dimension(4)    :: b
  real, dimension(4)    :: phi
  integer, dimension(4) :: vert
  integer               :: i, j, k, e, c, s
  real                  :: x, y, z, x2, y2, z2, xy, xz, yz, rho
!==============================================================================!

  !----------------------------------------------!
  !   Find curvatures at sides, and distribute   !
  !      them to elements surrounding them       !
  !----------------------------------------------!

  surf % elem(1:surf % n_elems) % curv = 0.0
  surf % elem(1:surf % n_elems) % xc   = 0.0
  surf % elem(1:surf % n_elems) % yc   = 0.0
  surf % elem(1:surf % n_elems) % zc   = 0.0

  do s = 1, surf % n_sides
    vert(1) = surf % side(s) % a
    vert(2) = surf % side(s) % b
    vert(3) = surf % side(s) % c
    vert(4) = surf % side(s) % d
    a(:,:) = 0
    b(:)   = 0

    do i = 1, 4

      x = surf % vert(vert(i)) % x_n
      y = surf % vert(vert(i)) % y_n
      z = surf % vert(vert(i)) % z_n

      x2 = x * x;  y2 = y * y;  z2 = z * z
      xy = x * y;  xz = x * z;  yz = y * z

      ! Fill the diagonal
      a(1,1) = a(1,1) + 1;  a(2,2) = a(2,2) + x2
      a(3,3) = a(3,3) + y2; a(4,4) = a(4,4) + z2

      ! Off-diagonal terms withouth first row and column
      a(2,3) = a(2,3) + xy;  a(2,4) = a(2,4) + xz;  a(3,4) = a(3,4) + yz;
      a(3,2) = a(3,2) + xy;  a(4,2) = a(4,2) + xz;  a(4,3) = a(4,3) + yz

      ! First row and column
      a(1,2) = a(1,2) + x;  a(1,3) = a(1,3) + y;  a(1,4) = a(1,4) + z;
      a(2,1) = a(2,1) + x;  a(3,1) = a(3,1) + y;  a(4,1) = a(4,1) + z

      ! Right hand side
      rho = x2 + y2 + z2
      b(1) = b(1) - rho
      b(2) = b(2) - rho * x;  b(3) = b(3) - rho * y;  b(4) = b(4) - rho * z
    end do

    call Math_Mod_Gaussian_Elimination(4, a, b, phi)

    ! Center of the sphere (could be stored in elems / verts too)
    x = -0.5 * phi(2);  y = -0.5 * phi(3);  z = -0.5 * phi(4)

    ! Sphere radius
    rho = sqrt(x*x + y*y + z*z - phi(1))

    surf % elem(surf % side(s) % ea) % curv =  &
    surf % elem(surf % side(s) % ea) % curv + 1.0 / rho * ONE_THIRD
    surf % elem(surf % side(s) % eb) % curv =  &
    surf % elem(surf % side(s) % eb) % curv + 1.0 / rho * ONE_THIRD

    ! This looks kind of crude
    surf % elem(surf % side(s) % ea) % xc =  &
    surf % elem(surf % side(s) % ea) % xc + x * ONE_THIRD
    surf % elem(surf % side(s) % ea) % yc =  &
    surf % elem(surf % side(s) % ea) % yc + y * ONE_THIRD
    surf % elem(surf % side(s) % ea) % zc =  &
    surf % elem(surf % side(s) % ea) % zc + z * ONE_THIRD
    surf % elem(surf % side(s) % eb) % xc =  &
    surf % elem(surf % side(s) % eb) % xc + x * ONE_THIRD
    surf % elem(surf % side(s) % eb) % yc =  &
    surf % elem(surf % side(s) % eb) % yc + y * ONE_THIRD
    surf % elem(surf % side(s) % eb) % zc =  &
    surf % elem(surf % side(s) % eb) % zc + z * ONE_THIRD
  end do

  ! Compute average curvature (for debugging)
  ! rho = 0
  ! do e = 1, surf % n_elems
  !   rho = rho + surf % elem(e) % curv
  ! end do
  ! rho = rho / surf % n_elems
  ! print *, 'average curvature = ', rho

  !-------------------------------------------------------------------!
  !   Interpolate normals at nodes from values in surrounding elems   !
  !-------------------------------------------------------------------!
  surf % vert(1:surf % n_verts) % curv = 0.
  do e = 1, surf % n_elems

    i = surf % elem(e) % i
    j = surf % elem(e) % j
    k = surf % elem(e) % k
    surf % vert(i) % curv = surf % vert(i) % curv  &
                          + surf % elem(e) % curv / real(surf % vert(i) % nne)
    surf % vert(j) % curv = surf % vert(j) % curv  &
                          + surf % elem(e) % curv / real(surf % vert(j) % nne)
    surf % vert(k) % curv = surf % vert(k) % curv  &
                          + surf % elem(e) % curv / real(surf % vert(k) % nne)
  end do

  end subroutine
