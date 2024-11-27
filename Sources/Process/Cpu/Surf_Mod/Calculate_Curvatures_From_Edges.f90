!==============================================================================!
  subroutine Calculate_Curvatures_From_Edges(Surf)
!------------------------------------------------------------------------------!
!>  The subroutine Calculate_Curvatures_From_Edges in the Surf_Mod module
!>  focuses on calculating surface curvatures from edges, involving four
!>  points, and subsequently distributing these values to the corresponding
!>  elements and nodes.  As the simplest and least accurate method among three
!>  available in the module, it employs a lower order of curvature calculation.
!>  This process involves constructing spheres at each side, determining their
!>  centers and radii, and averaging the curvature contributions to the
!>  elements and nodes.  Such curvature calculations are critical for
!>  estimating surface tension forces in VOF interface tracking methods.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf  !! parent class
!-----------------------------------[Locals]-----------------------------------!
  real, dimension(4,4)  :: a
  real, dimension(4)    :: b
  real, dimension(4)    :: phi
  integer, dimension(4) :: Vert
  integer               :: i, j, k, e, s
  logical               :: invertible
  real                  :: x, y, z, x2, y2, z2, xy, xz, yz, rho
!==============================================================================!

  !----------------------------------------------!
  !   Find curvatures at sides, and distribute   !
  !      them to elements surrounding them       !
  !----------------------------------------------!

  Surf % Elem(1:Surf % n_elems) % curv = 0.0
  Surf % Elem(1:Surf % n_elems) % xc   = 0.0
  Surf % Elem(1:Surf % n_elems) % yc   = 0.0
  Surf % Elem(1:Surf % n_elems) % zc   = 0.0

  do s = 1, Surf % n_sides
    Vert(1) = Surf % side(s) % a
    Vert(2) = Surf % side(s) % b
    Vert(3) = Surf % side(s) % c
    Vert(4) = Surf % side(s) % d
    a(:,:) = 0
    b(:)   = 0

    do i = 1, 4

      x = Surf % Vert(Vert(i)) % x_n
      y = Surf % Vert(Vert(i)) % y_n
      z = Surf % Vert(Vert(i)) % z_n

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

    call Math % Gaussian_Elimination(4, a, b, phi, invertible)

    ! If vertices are not co-planar
    if(invertible) then

      ! Center of the sphere (could be stored in elems / verts too)
      x = -0.5 * phi(2);  y = -0.5 * phi(3);  z = -0.5 * phi(4)

      ! Sphere radius
      rho = sqrt(x*x + y*y + z*z - phi(1))

    else
      rho = PETA  ! some big number
    end if

    Surf % Elem(Surf % side(s) % ea) % curv =  &
    Surf % Elem(Surf % side(s) % ea) % curv + 1.0 / rho * ONE_THIRD
    Surf % Elem(Surf % side(s) % eb) % curv =  &
    Surf % Elem(Surf % side(s) % eb) % curv + 1.0 / rho * ONE_THIRD

    ! This looks kind of crude
    Surf % Elem(Surf % side(s) % ea) % xc =  &
    Surf % Elem(Surf % side(s) % ea) % xc + x * ONE_THIRD
    Surf % Elem(Surf % side(s) % ea) % yc =  &
    Surf % Elem(Surf % side(s) % ea) % yc + y * ONE_THIRD
    Surf % Elem(Surf % side(s) % ea) % zc =  &
    Surf % Elem(Surf % side(s) % ea) % zc + z * ONE_THIRD
    Surf % Elem(Surf % side(s) % eb) % xc =  &
    Surf % Elem(Surf % side(s) % eb) % xc + x * ONE_THIRD
    Surf % Elem(Surf % side(s) % eb) % yc =  &
    Surf % Elem(Surf % side(s) % eb) % yc + y * ONE_THIRD
    Surf % Elem(Surf % side(s) % eb) % zc =  &
    Surf % Elem(Surf % side(s) % eb) % zc + z * ONE_THIRD
  end do

  ! Compute average curvature (for debugging)
  ! rho = 0
  ! do e = 1, Surf % n_elems
  !   rho = rho + Surf % Elem(e) % curv
  ! end do
  ! rho = rho / Surf % n_elems
  ! print *, 'average curvature = ', rho

  !-------------------------------------------------------------------!
  !   Interpolate normals at nodes from values in surrounding elems   !
  !-------------------------------------------------------------------!
  Surf % Vert(1:Surf % n_verts) % curv = 0.
  do e = 1, Surf % n_elems

    i = Surf % Elem(e) % v(1)
    j = Surf % Elem(e) % v(2)
    k = Surf % Elem(e) % v(3)
    Surf % Vert(i) % curv = Surf % Vert(i) % curv  &
                          + Surf % Elem(e) % curv / real(Surf % Vert(i) % nne)
    Surf % Vert(j) % curv = Surf % Vert(j) % curv  &
                          + Surf % Elem(e) % curv / real(Surf % Vert(j) % nne)
    Surf % Vert(k) % curv = Surf % Vert(k) % curv  &
                          + Surf % Elem(e) % curv / real(Surf % Vert(k) % nne)
  end do

  end subroutine
