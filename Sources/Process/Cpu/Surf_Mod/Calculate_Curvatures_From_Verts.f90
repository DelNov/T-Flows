!==============================================================================!
  subroutine Calculate_Curvatures_From_Verts(Surf)
!------------------------------------------------------------------------------!
!>  The subroutine Calculate_Curvatures_From_Verts in the Surf_Mod module
!>  calculates surface curvatures by considering all neighboring vertices
!>  around a given vertex, typically ranging from five to seven neighbors.
!>  This method represents an intermediate level of curvature calculation
!>  in terms of complexity and accuracy. It employs a more sophisticated
!>  approach than Calculate_Curvatures_From_Edges by using a larger set of
!>  points for curvature estimation. This increased detail allows for more
!>  accurate curvature values, which are crucial for correctly modeling
!>  surface tension in VOF interface tracking. The calculated curvatures
!>  are then distributed to the elements. This enhanced approach is
!>  significant for simulations where precise curvature estimations impact
!>  the physical accuracy of the modeled phenomena.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf  !! parent class
!-----------------------------------[Locals]-----------------------------------!
  real, dimension(4,4) :: a
  real, dimension(4)   :: b
  real, dimension(4)   :: phi
  integer, allocatable :: vert_v(:,:)
  integer              :: v, i, j, k, e, c, d, s, max_nnv
  logical              :: invertible
  real                 :: x, y, z, x2, y2, z2, xy, xz, yz, rho
!==============================================================================!

  ! Work out number of vertices around each vertex
  Surf % Vert(1:Surf % n_verts) % nnv = 0
  do s = 1, Surf % n_sides
    c = Surf % side(s) % c
    d = Surf % side(s) % d
    Surf % Vert(c) % nnv = Surf % Vert(c) % nnv + 1
    Surf % Vert(d) % nnv = Surf % Vert(d) % nnv + 1
  end do

  max_nnv = maxval(Surf % Vert(1:Surf % n_verts) % nnv)

  ! Form vert_v structure
  allocate( vert_v(max_nnv, Surf % n_verts) );  vert_v = 0

  Surf % Vert(1:Surf % n_verts) % nnv = 0
  do s = 1, Surf % n_sides
    c = Surf % side(s) % c
    d = Surf % side(s) % d

    Surf % Vert(c) % nnv = Surf % Vert(c) % nnv + 1
    vert_v(Surf % Vert(c) % nnv, c) = d

    Surf % Vert(d) % nnv = Surf % Vert(d) % nnv + 1
    vert_v(Surf % Vert(d) % nnv, d) = c
  end do

  !----------------------------------------------!
  !   Find curvatures at nodes, and distribute   !
  !      them to elements surrounding them       !
  !----------------------------------------------!

  Surf % Elem(1:Surf % n_elems) % curv = 0.0
  Surf % Elem(1:Surf % n_elems) % xc   = 0.0
  Surf % Elem(1:Surf % n_elems) % yc   = 0.0
  Surf % Elem(1:Surf % n_elems) % zc   = 0.0

  do v = 1, Surf % n_verts

    a(:,:) = 0
    b(:)   = 0

    do k = 1, Surf % Vert(v) % nnv

      x = Surf % Vert( vert_v(k,v) ) % x_n
      y = Surf % Vert( vert_v(k,v) ) % y_n
      z = Surf % Vert( vert_v(k,v) ) % z_n

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

    Surf % Vert(v) % curv = 1.0 / rho
  end do

  ! Compute average curvature (for debugging)
  ! rho = 0
  ! do v = 1, Surf % n_verts
  !   rho = rho + Surf % Vert(v) % curv
  ! end do
  ! rho = rho / Surf % n_elems
  ! print *, 'average curvature = ', rho

  !----------------------------------------------------------------------!
  !   Interpolate normals at elems from values in surrounding vertices   !
  !----------------------------------------------------------------------!
  Surf % Elem(1:Surf % n_elems) % curv = 0.
  do e = 1, Surf % n_elems

    i = Surf % Elem(e) % v(1)
    j = Surf % Elem(e) % v(2)
    k = Surf % Elem(e) % v(3)
    Surf % Elem(e) % curv = ONE_THIRD * (  Surf % Vert(i) % curv  &
                                         + Surf % Vert(j) % curv  &
                                         + Surf % Vert(k) % curv )
  end do

  end subroutine
