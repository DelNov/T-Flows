!==============================================================================!
  subroutine Front_Mod_Calculate_Curvatures_From_Elems(front)
!------------------------------------------------------------------------------!
!   Calculates surface curvatures from elements (using its own vertices and    !
!   the vertices from neighbouring elements, around twelve points) and         !
!   distributes the valies to the vertices.                                    !
!                                                                              !
!   This is the highest order of curvature calculation in the code.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Front_Type), target :: front
!-----------------------------------[Locals]-----------------------------------!
  real, dimension(4,4)  :: a
  real, dimension(4)    :: b
  real, dimension(4)    :: phi
  integer, dimension(4) :: node
  integer, allocatable  :: elem_n_verts(:)
  integer, allocatable  :: vert_v(:,:)
  integer, allocatable  :: elem_v(:,:)
  integer               :: v, i, j, k, e, c, d, s, max_nnv
  real                  :: x, y, z, x2, y2, z2, xy, xz, yz, rho
!==============================================================================!

  ! Work out number of vertices around each vertex
  front % vert(1:front % n_verts) % nnv = 0
  do s = 1, front % n_sides
    c = front % side(s) % c
    d = front % side(s) % d
    front % vert(c) % nnv = front % vert(c) % nnv + 1
    front % vert(d) % nnv = front % vert(d) % nnv + 1
  end do

  max_nnv = maxval(front % vert(1:front % n_verts) % nnv)

  ! Form vert_v structure
  allocate( vert_v(max_nnv, front % n_verts) );  vert_v = 0

  front % vert(1:front % n_verts) % nnv = 0
  do s = 1, front % n_sides
    c = front % side(s) % c
    d = front % side(s) % d

    front % vert(c) % nnv = front % vert(c) % nnv + 1
    vert_v(front % vert(c) % nnv, c) = d

    front % vert(d) % nnv = front % vert(d) % nnv + 1
    vert_v(front % vert(d) % nnv, d) = c
  end do

  ! Work out (maximum) number of vertices around each element
  allocate(elem_n_verts(front % n_elems))
  elem_n_verts(:) = 0
  do e = 1, front % n_elems
    i = front % elem(e) % i
    j = front % elem(e) % j
    k = front % elem(e) % k
    elem_n_verts(e) = elem_n_verts(e) + front % vert(i) % nnv
    elem_n_verts(e) = elem_n_verts(e) + front % vert(j) % nnv
    elem_n_verts(e) = elem_n_verts(e) + front % vert(k) % nnv
  end do

  max_nnv = maxval(elem_n_verts(1:front % n_elems))
  allocate(elem_v(max_nnv, front % n_elems))
  elem_v(:,:) = HUGE_INT

  elem_n_verts(:) = 0
  do e = 1, front % n_elems
    i = front % elem(e) % i
    j = front % elem(e) % j
    k = front % elem(e) % k
    do v = 1, front % vert(i) % nnv
      elem_n_verts(e) = elem_n_verts(e) + 1;
      elem_v(elem_n_verts(e), e) = vert_v(v, i)
    end do
    do v = 1, front % vert(j) % nnv
      elem_n_verts(e) = elem_n_verts(e) + 1;
      elem_v(elem_n_verts(e), e) = vert_v(v, j)
    end do
    do v = 1, front % vert(k) % nnv
      elem_n_verts(e) = elem_n_verts(e) + 1;
      elem_v(elem_n_verts(e), e) = vert_v(v, k)
    end do
  end do

  ! Remove duplicate entries in the vertex list
  do e = 1, front % n_elems
    call Sort_Mod_Unique_Int(elem_v(1:elem_n_verts(e), e), elem_n_verts(e))
  end do

  !---------------------------------------------!
  !   Find curvatures at elements using first   !
  !        and second neighbouring nodes        !
  !---------------------------------------------!

  front % elem(1:front % n_elems) % curv = 0.0
  front % elem(1:front % n_elems) % xc   = 0.0
  front % elem(1:front % n_elems) % yc   = 0.0
  front % elem(1:front % n_elems) % zc   = 0.0

  do e = 1, front % n_elems

    a(:,:) = 0
    b(:)   = 0

    do k = 1, elem_n_verts(e)

      x = front % vert( elem_v(k,e) ) % x_n
      y = front % vert( elem_v(k,e) ) % y_n
      z = front % vert( elem_v(k,e) ) % z_n

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

    front % elem(e) % curv = 1.0 / rho
    front % elem(e) % xc   = x
    front % elem(e) % yc   = y
    front % elem(e) % zc   = z
  end do

  ! Compute average curvature (for debugging)
  ! rho = 0
  ! do e = 1, front % n_elems
  !   rho = rho + front % elem(e) % curv
  ! end do
  ! rho = rho / front % n_elems
  ! print *, 'average curvature = ', rho

  !-------------------------------------------------------------------!
  !   Interpolate normals at nodes from values in surrounding elems   !
  !-------------------------------------------------------------------!
  front % vert(1:front % n_verts) % curv = 0.
  do e = 1, front % n_elems

    i = front % elem(e) % i
    j = front % elem(e) % j
    k = front % elem(e) % k
    front % vert(i) % curv = front % vert(i) % curv  &
                           + front % elem(e) % curv/real(front % vert(i) % nne)
    front % vert(j) % curv = front % vert(j) % curv  &
                           + front % elem(e) % curv/real(front % vert(j) % nne)
    front % vert(k) % curv = front % vert(k) % curv  &
                           + front % elem(e) % curv/real(front % vert(k) % nne)
  end do

  end subroutine
