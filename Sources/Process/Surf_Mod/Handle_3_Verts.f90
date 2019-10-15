!==============================================================================!
  subroutine Surf_Mod_Handle_3_Verts(surf, surf_v)
!------------------------------------------------------------------------------!
!   Places surface where variable phi has value val                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Surf_Type), target :: surf
  real, dimension(3)      :: surf_v
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: vert(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ne
  integer                  :: ver(3)
  real                     :: a(3), b(3), tri_v(3)
!==============================================================================!

  ! Take aliases
  nv   => surf % n_verts
  ne   => surf % n_elems
  vert => surf % vert
  elem => surf % elem

  ne = ne + 1           ! one new element

  ! Store last three vertices for the new element
  elem(ne) % i = nv - 2
  elem(ne) % j = nv - 1
  elem(ne) % k = nv

  ! Take vectors a and b; u connection points 2 and 1;
  ! v connecting 3 and 1
  ver(1) = elem(ne) % i
  ver(2) = elem(ne) % j
  ver(3) = elem(ne) % k

  a(1) = vert(ver(2)) % x_n - vert(ver(1)) % x_n
  a(2) = vert(ver(2)) % y_n - vert(ver(1)) % y_n
  a(3) = vert(ver(2)) % z_n - vert(ver(1)) % z_n

  b(1) = vert(ver(3)) % x_n - vert(ver(1)) % x_n
  b(2) = vert(ver(3)) % y_n - vert(ver(1)) % y_n
  b(3) = vert(ver(3)) % z_n - vert(ver(1)) % z_n

  tri_v = Math_Mod_Cross_Product(a, b)

  if(dot_product(surf_v, tri_v) < 0.0) then
    call Swap_Int(elem(ne) % j,  &
                  elem(ne) % k)
  end if

  end subroutine
