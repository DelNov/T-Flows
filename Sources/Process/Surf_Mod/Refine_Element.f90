!==============================================================================!
  subroutine Surf_Mod_Refine_Element(surf, e)
!------------------------------------------------------------------------------!
!   Create new vertex in the center of existing element specified with e       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Surf_Type), target :: surf
  integer                 :: e
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: vert(:)
  type(Side_Type), pointer :: side(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ns, ne
  integer                  :: i, j, k, ei, ej, ek, si, sj, sk
!==============================================================================!

  ! Take aliases
  nv   => surf % n_verts
  ns   => surf % n_sides
  ne   => surf % n_elems
  vert => surf % vert
  side => surf % side
  elem => surf % elem

  i = elem(e) % i
  j = elem(e) % j
  k = elem(e) % k

  ei = elem(e) % ei
  ej = elem(e) % ej
  ek = elem(e) % ek

  si = elem(e) % si
  sj = elem(e) % sj
  sk = elem(e) % sk

  !-----------------------------------------------!
  !   Create new vertex in the element's center   !
  !-----------------------------------------------!
  nv = nv+ 1
  vert(nv) % x_n = (vert(i) % x_n + vert(j) % x_n + vert(k) % x_n) * ONE_THIRD
  vert(nv) % y_n = (vert(i) % y_n + vert(j) % y_n + vert(k) % y_n) * ONE_THIRD
  vert(nv) % z_n = (vert(i) % z_n + vert(j) % z_n + vert(k) % z_n) * ONE_THIRD

  !----------------------------!
  !   Create three new sides   !
  !----------------------------!
  ns = ns + 1
  side(ns) % c = i
  side(ns) % d = nv
  side(ns) % a = k
  side(ns) % b = j
  side(ns) % ea = ne + 2  !  yet to be formed
  side(ns) % eb = e

  ns = ns + 1
  side(ns) % c = j
  side(ns) % d = nv
  side(ns) % a = i
  side(ns) % b = k
  side(ns) % ea = e
  side(ns) % eb = ne + 1  ! yet to be formed

  ns = ns + 1
  side(ns) % c = k
  side(ns) % d = nv
  side(ns) % a = j
  side(ns) % b = i
  side(ns) % ea = ne + 1  ! yet to be formed
  side(ns) % eb = ne + 2  ! yet to be formed

  !-------------------------------!
  !   Create three new elements   !
  !-------------------------------!

  ! Handle new elements' nodes
  elem(e) % k = nv   ! this element actually existed

  ne = ne + 1
  elem(ne) % i = nv
  elem(ne) % j = j
  elem(ne) % k = k

  ne = ne + 1
  elem(ne) % i = i
  elem(ne) % j = nv
  elem(ne) % k = k

  ! Handle new elements' neighbours
  elem(e) % ei = ne - 1
  elem(e) % ek = ne

  elem(ne-1) % ei = ei
  elem(ne-1) % ej = ej
  elem(ne-1) % ek = e

  elem(ne) % ei = ne - 1
  elem(ne) % ej = ej
  elem(ne) % ek = ek

  ! Handle new elements' sides
  elem(e) % si = ns - 1     ! old j and new
  elem(e) % sj = ns - 2     ! old i and new
  elem(e) % sk = sk

  elem(ne-1) % si = si
  elem(ne-1) % sj = ns
  elem(ne-1) % sk = ns - 1  ! old j and new

  elem(ne) % si = ns        ! old k and new
  elem(ne) % sj = sj
  elem(ne) % sk = ns - 2    ! old i and new

  end subroutine
