!==============================================================================!
  subroutine Surf_Mod_Calculate_Curvatures_From_Verts(surf)
!------------------------------------------------------------------------------!
!   Calculates surface curvatures from vertices (all neighbous, five to seven  !
!   on the average) and distributes the values to elements.                    !
!                                                                              !
!   This is the intermediate order of curvature calculation in the code.       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Surf_Type), target :: surf
!-----------------------------------[Locals]-----------------------------------!
  real, dimension(4,4)  :: a
  real, dimension(4)    :: b
  real, dimension(4)    :: phi
  integer, dimension(4) :: node
  integer, allocatable  :: vert_v(:,:)
  integer               :: v, i, j, k, e, c, d, s, max_nnv
  real                  :: x, y, z, x2, y2, z2, xy, xz, yz, rho
!==============================================================================!

  ! Work out number of vertices around each vertex
  surf % vert(1:surf % n_verts) % nnv = 0
  do s = 1, surf % n_sides
    c = surf % side(s) % c
    d = surf % side(s) % d
    surf % vert(c) % nnv = surf % vert(c) % nnv + 1
    surf % vert(d) % nnv = surf % vert(d) % nnv + 1
  end do

  max_nnv = maxval(surf % vert(1:surf % n_verts) % nnv)

  ! Form vert_v structure
  allocate( vert_v(max_nnv, surf % n_verts) );  vert_v = 0

  surf % vert(1:surf % n_verts) % nnv = 0
  do s = 1, surf % n_sides
    c = surf % side(s) % c
    d = surf % side(s) % d

    surf % vert(c) % nnv = surf % vert(c) % nnv + 1
    vert_v(surf % vert(c) % nnv, c) = d

    surf % vert(d) % nnv = surf % vert(d) % nnv + 1
    vert_v(surf % vert(d) % nnv, d) = c
  end do

  !----------------------------------------------!
  !   Find curvatures at nodes, and distribute   !
  !      them to elements surrounding them       !
  !----------------------------------------------!

  surf % elem(1:surf % n_elems) % curv = 0.0
  surf % elem(1:surf % n_elems) % xc   = 0.0
  surf % elem(1:surf % n_elems) % yc   = 0.0
  surf % elem(1:surf % n_elems) % zc   = 0.0

  do v = 1, surf % n_verts

    a(:,:) = 0
    b(:)   = 0

    do k = 1, surf % vert(v) % nnv

      x = surf % vert( vert_v(k,v) ) % x_n
      y = surf % vert( vert_v(k,v) ) % y_n
      z = surf % vert( vert_v(k,v) ) % z_n

      x2 = x * x;  y2 = y * y;  z2 = z * z
      xy = x * y;  xz = x * z;  yz = y * z

      ! Fill the diagonal
      a(1,1) = a(1,1) + 1
      a(2,2) = a(2,2) + x2;  a(3,3) = a(3,3) + y2; a(4,4) = a(4,4) + z2

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

    surf % vert(v) % curv = 1.0 / rho
  end do

  ! Compute average curvature (for debugging)
  ! rho = 0
  ! do v = 1, surf % n_verts
  !   rho = rho + surf % vert(v) % curv
  ! end do
  ! rho = rho / surf % n_elems
  ! print *, 'average curvature = ', rho

  !----------------------------------------------------------------------!
  !   Interpolate normals at elems from values in surrounding vertices   !
  !----------------------------------------------------------------------!
  surf % elem(1:surf % n_elems) % curv = 0.
  do e = 1, surf % n_elems

    i = surf % elem(e) % v(1)
    j = surf % elem(e) % v(2)
    k = surf % elem(e) % v(3)
    surf % elem(e) % curv = ONE_THIRD * (  surf % vert(i) % curv  &
                                         + surf % vert(j) % curv  &
                                         + surf % vert(k) % curv )
  end do

  end subroutine
