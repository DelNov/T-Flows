!==============================================================================!
  subroutine Handle_5_Points(Surf, surf_v, enforce_triangles)
!------------------------------------------------------------------------------!
!>  The Handle_5_Points subroutine efficiently manages the creation of surface
!>  elements when an interface intersects a cell at five points. It determines
!>  the optimal arrangement of these points to form either a single
!>  pentagonal element or three triangular elements, based on the specified
!>  enforcement of triangular elements.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Alias and variable setup: Establishes aliases for the number of          !
!     vertices (nv), elements (ne), and pointers to vertices (Vert) and        !
!     elements (Elem). Also, sets up permutations for vertex orderings.        !
!   * Permutation testing:                                                     !
!     - Iterates through different permutations of the five intersecting       !
!       vertices to find an arrangement where the vertices form a coherent     !
!       surface element(s).                                                    !
!     - For each permutation, calculates vectors between vertices and their    !
!       cross products to determine the surface normals.                       !
!   * Orientation checking:                                                    !
!     - Checks the orientation of the potential elements by comparing the dot  !
!       product of the surface vector (surf_v) with the calculated normals.    !
!     - Selects the permutation where the normals align correctly with the     !
!       surface vector.                                                        !
!   * Element formation:                                                       !
!     - Depending on the enforce_triangles flag:                               !
!       > If not enforcing triangles, creates one pentagonal element with      !
!         the five vertices.                                                   !
!       > If enforcing triangles, creates three triangular elements from the   !
!         five vertices.                                                       !
!     - Assigns vertices to the new element(s) based on the successful         !
!       permutation.                                                           !
!   * Error Handling:                                                          !
!     - Includes a message and error handling if no suitable permutation is    !
!       found, indicating a critical issue in the surface generation process.  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf               !! parent class
  real                     :: surf_v(3)          !! surface vector
  logical                  :: enforce_triangles  !! controls if creation of
                                                 !! triangles is enforced
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: Vert(:)
  type(Elem_Type), pointer :: Elem(:)
  integer,         pointer :: nv, ne
  integer                  :: ver(5), loop
  real                     :: v_21(3), v_31(3), v_41(3), v_51(3)
  real                     :: tri_p_123(3), tri_p_134(3), tri_p_145(3)
  integer                  :: permutations(5, 24)
!=============================================================================!

  ! Take aliases
  nv   => Surf % n_verts
  ne   => Surf % n_elems
  Vert => Surf % Vert
  Elem => Surf % Elem

  permutations = reshape((/ &
    0, 1, 2, 3, 4,  &
    0, 2, 1, 3, 4,  &
    0, 3, 1, 2, 4,  &
    0, 1, 3, 2, 4,  &
    0, 2, 3, 1, 4,  &
    0, 3, 2, 1, 4,  &
    0, 4, 2, 1, 3,  &
    0, 2, 4, 1, 3,  &
    0, 1, 4, 2, 3,  &
    0, 4, 1, 2, 3,  &
    0, 2, 1, 4, 3,  &
    0, 1, 2, 4, 3,  &
    0, 1, 3, 4, 2,  &
    0, 3, 1, 4, 2,  &
    0, 4, 1, 3, 2,  &
    0, 1, 4, 3, 2,  &
    0, 3, 4, 1, 2,  &
    0, 4, 3, 1, 2,  &
    0, 4, 3, 2, 1,  &
    0, 3, 4, 2, 1,  &
    0, 2, 4, 3, 1,  &
    0, 4, 2, 3, 1,  &
    0, 3, 2, 4, 1,  &
    0, 2, 3, 4, 1   &
  /), shape(permutations))

  !-------------------------------------------!
  !   Try to find a permutation which works   !
  !-------------------------------------------!
  do loop = 1, 24

    ver(1) = nv - permutations(1, loop)
    ver(2) = nv - permutations(2, loop)
    ver(3) = nv - permutations(3, loop)
    ver(4) = nv - permutations(4, loop)
    ver(5) = nv - permutations(5, loop)

    v_21(1) = Vert(ver(2)) % x_n - Vert(ver(1)) % x_n
    v_21(2) = Vert(ver(2)) % y_n - Vert(ver(1)) % y_n
    v_21(3) = Vert(ver(2)) % z_n - Vert(ver(1)) % z_n

    v_31(1) = Vert(ver(3)) % x_n - Vert(ver(1)) % x_n
    v_31(2) = Vert(ver(3)) % y_n - Vert(ver(1)) % y_n
    v_31(3) = Vert(ver(3)) % z_n - Vert(ver(1)) % z_n

    v_41(1) = Vert(ver(4)) % x_n - Vert(ver(1)) % x_n
    v_41(2) = Vert(ver(4)) % y_n - Vert(ver(1)) % y_n
    v_41(3) = Vert(ver(4)) % z_n - Vert(ver(1)) % z_n

    v_51(1) = Vert(ver(5)) % x_n - Vert(ver(1)) % x_n
    v_51(2) = Vert(ver(5)) % y_n - Vert(ver(1)) % y_n
    v_51(3) = Vert(ver(5)) % z_n - Vert(ver(1)) % z_n

    tri_p_123 = Math % Cross_Product(v_21, v_31)
    tri_p_134 = Math % Cross_Product(v_31, v_41)
    tri_p_145 = Math % Cross_Product(v_41, v_51)

    if(dot_product(surf_v, tri_p_123) > 0.0 .and.  &
       dot_product(surf_v, tri_p_134) > 0.0 .and.  &
       dot_product(surf_v, tri_p_145) > 0.0) then
      goto 1
    end if

  end do

  call Message % Error(72,                                                    &
                  'Failed to find a good permutation in Handle_5_Points! '//  &
                  'This error is critical.  Exiting!',                        &
                  file=__FILE__, line=__LINE__)

  !--------------------------------------------------------------------!
  !   You've found a permutation which works, store new elements now   !
  !--------------------------------------------------------------------!
1 continue

  if(.not. enforce_triangles) then

    ! One new element with five vertices
    ne = ne + 1
    Elem(ne) % nv = 5
    Elem(ne) % v(1:5) = ver(1:5)

  else  ! only triangles are allowed

    ! One new element with three vertices
    ne = ne + 1
    Elem(ne) % nv = 3
    Elem(ne) % v(1) = ver(1)
    Elem(ne) % v(2) = ver(2)
    Elem(ne) % v(3) = ver(3)

    ! Second new element with three vertices
    ne = ne + 1
    Elem(ne) % nv = 3
    Elem(ne) % v(1) = ver(1)
    Elem(ne) % v(2) = ver(3)
    Elem(ne) % v(3) = ver(4)

    ! Third new element with three vertices
    ne = ne + 1
    Elem(ne) % nv = 3
    Elem(ne) % v(1) = ver(1)
    Elem(ne) % v(2) = ver(4)
    Elem(ne) % v(3) = ver(5)

  end if

  end subroutine