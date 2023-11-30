!==============================================================================!
  subroutine Calculate_Curvatures_From_Elems(Surf)
!------------------------------------------------------------------------------!
!>  The subroutine Calculate_Curvatures_From_Elems in the Surf_Mod module is
!>  the most advanced and accurate method for calculating surface curvatures
!>  in the code. It computes curvatures using each element's vertices and
!>  those of neighboring elements, typically involving around twelve points.
!>  This method offers the highest level of curvature calculation accuracy by
!>  utilizing a larger, more comprehensive set of points. These calculated
!>  curvatures are then distributed to the vertices. The precision of this
!>  method is crucial for accurately modeling surface phenomena, particularly
!>  in simulations where surface tension forces are critical.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf  !! parent class
!-----------------------------------[Locals]-----------------------------------!
  real, dimension(4,4) :: a
  real, dimension(4)   :: b
  real, dimension(4)   :: phi
  integer, allocatable :: elem_n_verts(:)
  integer, allocatable :: vert_v(:,:)
  integer, allocatable :: elem_v(:,:)
  integer              :: v, i_ver, j_ver, e, c, d, s, max_nnv
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

  ! Work out (maximum) number of vertices around each element
  allocate(elem_n_verts(Surf % n_elems))
  elem_n_verts(:) = 0
  do e = 1, Surf % n_elems
    do i_ver = 1, Surf % Elem(e) % nv
      v = Surf % Elem(e) % v(i_ver)
      elem_n_verts(e) = elem_n_verts(e) + Surf % Vert(v) % nnv
    end do
  end do

  max_nnv = maxval(elem_n_verts(1:Surf % n_elems))
  allocate(elem_v(max_nnv, Surf % n_elems))
  elem_v(:,:) = HUGE_INT

  elem_n_verts(:) = 0
  do e = 1, Surf % n_elems
    do i_ver = 1, Surf % Elem(e) % nv
      v = Surf % Elem(e) % v(i_ver)
      do j_ver = 1, Surf % Vert(v) % nnv
        elem_n_verts(e) = elem_n_verts(e) + 1;
        elem_v(elem_n_verts(e), e) = vert_v(j_ver, v)
      end do
    end do
  end do

  ! Remove duplicate entries in the vertex list
  do e = 1, Surf % n_elems
    call Sort % Unique_Int(elem_v(1:elem_n_verts(e), e), elem_n_verts(e))
  end do

  !---------------------------------------------!
  !   Find curvatures at elements using first   !
  !        and second neighbouring nodes        !
  !---------------------------------------------!

  Surf % Elem(1:Surf % n_elems) % curv = 0.0
  Surf % Elem(1:Surf % n_elems) % xc   = 0.0
  Surf % Elem(1:Surf % n_elems) % yc   = 0.0
  Surf % Elem(1:Surf % n_elems) % zc   = 0.0

  do e = 1, Surf % n_elems

    a(:,:) = 0
    b(:)   = 0

    do i_ver = 1, elem_n_verts(e)

      x = Surf % Vert( elem_v(i_ver, e) ) % x_n
      y = Surf % Vert( elem_v(i_ver, e) ) % y_n
      z = Surf % Vert( elem_v(i_ver, e) ) % z_n

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
      x = 0.0;  y = 0.0;  z = 0.0  ! this is not super-smart, but OK
      rho = PETA                   ! some big number
    end if

    Surf % Elem(e) % curv = 1.0 / rho
    Surf % Elem(e) % xc   = x
    Surf % Elem(e) % yc   = y
    Surf % Elem(e) % zc   = z
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

    do i_ver = 1, Surf % Elem(e) % nv
      v = Surf % Elem(e) % v(i_ver)
      Surf % Vert(v) % curv = Surf % Vert(v) % curv  &
                      + Surf % Elem(e) % curv/real(Surf % Vert(v) % nne)
    end do
  end do

  end subroutine
