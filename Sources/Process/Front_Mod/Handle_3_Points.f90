!==============================================================================!
  subroutine Handle_3_Points(Front, surf_v, enforce_triangles)
!------------------------------------------------------------------------------!
!   Surface intersects cell at three points                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type), target :: Front
  real, dimension(3)        :: surf_v
  logical                   :: enforce_triangles
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: Vert(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ne
  integer                  :: ver(3)
  real                     :: a(3), b(3), tri_p(3)
!==============================================================================!

  ! Take aliases
  nv   => Front % n_verts
  ne   => Front % n_elems
  Vert => Front % Vert
  elem => Front % elem

  ! One new element with three vertices
  ne = ne + 1
  elem(ne) % nv = 3

  ! Store last three vertices for the new element
  elem(ne) % v(1) = nv - 2
  elem(ne) % v(2) = nv - 1
  elem(ne) % v(3) = nv

  ! Take vectors a and b; u connection points 2 and 1;
  ! v connecting 3 and 1
  ver(1) = elem(ne) % v(1)
  ver(2) = elem(ne) % v(2)
  ver(3) = elem(ne) % v(3)

  a(1) = Vert(ver(2)) % x_n - Vert(ver(1)) % x_n
  a(2) = Vert(ver(2)) % y_n - Vert(ver(1)) % y_n
  a(3) = Vert(ver(2)) % z_n - Vert(ver(1)) % z_n

  b(1) = Vert(ver(3)) % x_n - Vert(ver(1)) % x_n
  b(2) = Vert(ver(3)) % y_n - Vert(ver(1)) % y_n
  b(3) = Vert(ver(3)) % z_n - Vert(ver(1)) % z_n

  tri_p = Math_Mod_Cross_Product(a, b)

  if(dot_product(surf_v, tri_p) < 0.0) then
    call Swap_Int(elem(ne) % v(2),  &
                  elem(ne) % v(3))
  end if

  end subroutine
