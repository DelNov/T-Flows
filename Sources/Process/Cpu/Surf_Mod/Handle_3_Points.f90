!==============================================================================!
  subroutine Handle_3_Points(Surf, surf_v)
!------------------------------------------------------------------------------!
!>  The Handle_3_Points subroutine in the Surf_Type is a specialized function
!>  designed to handle the scenario where a surface intersects a cell at three
!>  points, forming a triangular element.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Alias assignment: Sets up aliases for the number of vertices (nv), the   !
!     number of elements (ne), vertices (Vert), and elements (Elem) for easier !
!     reference and manipulation.                                              !
!   * Element creation:                                                        !
!     - Increments the element count and defines a new element with three      !
!       vertices.                                                              !
!     - Assigns the last three vertices in the surface mesh to the new element.!
!   * Vector calculation:                                                      !
!     - Computes vectors a and b based on the positions of the vertices.       !
!       Vector a connects points 2 and 1; vector b connects points 3 and 1.    !
!     - Calculates the cross product of vectors a and b to determine the       !
!       orientation of the triangle in space.                                  !
!   * Orientation adjustment:                                                  !
!     - Checks the orientation of the newly formed triangle relative to the    !
!       surface vector (surf_v).                                               !
!     - If the dot product of the surface vector and the triangle's normal     !
!       (from the cross product) is negative, the subroutine adjusts the       !
!       vertex order to ensure correct orientation.                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf    !! parent class
  real, dimension(3)       :: surf_v  !! surface vector
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: Vert(:)
  type(Elem_Type), pointer :: Elem(:)
  integer,         pointer :: nv, ne
  integer                  :: ver(3)
  real                     :: a(3), b(3), tri_p(3)
!==============================================================================!

  ! Take aliases
  nv   => Surf % n_verts
  ne   => Surf % n_elems
  Vert => Surf % Vert
  Elem => Surf % Elem

  ! One new element with three vertices
  ne = ne + 1
  Elem(ne) % nv = 3

  ! Store last three vertices for the new element
  Elem(ne) % v(1) = nv - 2
  Elem(ne) % v(2) = nv - 1
  Elem(ne) % v(3) = nv

  ! Take vectors a and b; u connection points 2 and 1;
  ! v connecting 3 and 1
  ver(1) = Elem(ne) % v(1)
  ver(2) = Elem(ne) % v(2)
  ver(3) = Elem(ne) % v(3)

  a(1) = Vert(ver(2)) % x_n - Vert(ver(1)) % x_n
  a(2) = Vert(ver(2)) % y_n - Vert(ver(1)) % y_n
  a(3) = Vert(ver(2)) % z_n - Vert(ver(1)) % z_n

  b(1) = Vert(ver(3)) % x_n - Vert(ver(1)) % x_n
  b(2) = Vert(ver(3)) % y_n - Vert(ver(1)) % y_n
  b(3) = Vert(ver(3)) % z_n - Vert(ver(1)) % z_n

  tri_p = Math % Cross_Product(a, b)

  if(dot_product(surf_v, tri_p) < 0.0) then
    call Swap_Int(Elem(ne) % v(2),  &
                  Elem(ne) % v(3))
  end if

  end subroutine
