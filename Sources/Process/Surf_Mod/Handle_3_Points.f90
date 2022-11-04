!==============================================================================!
  subroutine Handle_3_Points(Surf, surf_v)
!------------------------------------------------------------------------------!
!   Surface intersects cell at three points                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf
  real, dimension(3)       :: surf_v
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
