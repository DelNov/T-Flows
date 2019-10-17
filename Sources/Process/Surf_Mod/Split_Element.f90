!==============================================================================!
  subroutine Surf_Mod_Split_Element(surf, e)
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
  integer                  :: s_in, s_jn, s_kn, e_ijn, e_njk, e_ink
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
  nv = nv + 1
  vert(nv) % x_n = (vert(i) % x_n + vert(j) % x_n + vert(k) % x_n) * ONE_THIRD
  vert(nv) % y_n = (vert(i) % y_n + vert(j) % y_n + vert(k) % y_n) * ONE_THIRD
  vert(nv) % z_n = (vert(i) % z_n + vert(j) % z_n + vert(k) % z_n) * ONE_THIRD

  ! Define a few handy constants
  s_in = ns + 1
  s_jn = ns + 2
  s_kn = ns + 3

  e_ijn = e
  e_njk = ne + 1
  e_ink = ne + 2

  !----------------------------!
  !   Create three new sides   !
  !----------------------------!
  side(s_in) % c  = i
  side(s_in) % d  = nv
  side(s_in) % a  = k
  side(s_in) % b  = j
  side(s_in) % ea = e_ink
  side(s_in) % eb = e_ijn

  side(s_jn) % c  = j
  side(s_jn) % d  = nv
  side(s_jn) % a  = i
  side(s_jn) % b  = k
  side(s_jn) % ea = e_ijn
  side(s_jn) % eb = e_njk

  side(s_kn) % c  = k
  side(s_kn) % d  = nv
  side(s_kn) % a  = j
  side(s_kn) % b  = i
  side(s_kn) % ea = e_njk
  side(s_kn) % eb = e_ink

  !-------------------------------!
  !   Create three new elements   !
  !-------------------------------!

  ! Handle new elements' nodes
  elem(e_ijn) % i = i
  elem(e_ijn) % j = j
  elem(e_ijn) % k = nv

  elem(e_njk) % i = nv
  elem(e_njk) % j = j
  elem(e_njk) % k = k

  elem(e_ink) % i = i
  elem(e_ink) % j = nv
  elem(e_ink) % k = k

  ! Handle new elements' neighbours
  elem(e_ijn) % ei = e_njk
  elem(e_ijn) % ej = e_ink
  elem(e_ijn) % ek = ek

  elem(e_njk) % ei = ei
  elem(e_njk) % ej = e_ink
  elem(e_njk) % ek = e_ijn

  elem(e_ink) % ei = e_njk
  elem(e_ink) % ej = ej
  elem(e_ink) % ek = e_ijn

  ! Handle new elements' sides
  elem(e_ijn) % si = s_jn
  elem(e_ijn) % sj = s_in
  elem(e_ijn) % sk = sk

  elem(e_njk) % si = si
  elem(e_njk) % sj = s_kn
  elem(e_njk) % sk = s_jn

  elem(e_ink) % si = s_kn
  elem(e_ink) % sj = sj
  elem(e_ink) % sk = s_in

  !-------------------------------------------!
  !   Increase number of sides and elements   !
  !-------------------------------------------!
  ns = ns + 3
  ne = ne + 2

  !-------------------------------------!
  !   Correct neighbours where needed   !
  !-------------------------------------!

  if(ei .gt. 0) then
    if(elem(ei) % ei .eq. e_ijn) elem(ei) % ei = e_njk
    if(elem(ei) % ej .eq. e_ijn) elem(ei) % ej = e_njk
    if(elem(ei) % ek .eq. e_ijn) elem(ei) % ek = e_njk
  end if

  if(side(si) % a .eq. i) side(si) % a = nv
  if(side(si) % b .eq. i) side(si) % b = nv

  if(ej .gt. 0) then
    if(elem(ej) % ei .eq. e_ijn) elem(ej) % ei = e_ink
    if(elem(ej) % ej .eq. e_ijn) elem(ej) % ej = e_ink
    if(elem(ej) % ek .eq. e_ijn) elem(ej) % ek = e_ink
  end if

  if(side(sj) % a .eq. j) side(si) % a = nv
  if(side(sj) % b .eq. j) side(si) % b = nv

  if(side(sk) % a .eq. k) side(sk) % a = nv
  if(side(sk) % b .eq. k) side(sk) % b = nv

!  vert(i)  % nne = vert(i) % nne + 1
!  vert(j)  % nne = vert(j) % nne + 1
!  vert(k)  % nne = vert(k) % nne + 1
!  vert(nv) % nne = 3

  end subroutine
