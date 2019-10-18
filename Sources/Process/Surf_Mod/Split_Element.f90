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
  integer                  :: i, j, k, ei, ej, ek, s, si, sj, sk
  integer                  :: sum_ijk, sum_cd, run, ea, eb
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

  if(ej .gt. 0) then
    if(elem(ej) % ei .eq. e_ijn) elem(ej) % ei = e_ink
    if(elem(ej) % ej .eq. e_ijn) elem(ej) % ej = e_ink
    if(elem(ej) % ek .eq. e_ijn) elem(ej) % ek = e_ink
  end if

  if(side(si) % ea .eq. e_ijn) then
    side(si) % ea = e_njk
    side(si) % a  = nv
  end if
  if(side(si) % eb .eq. e_ijn) then
    side(si) % eb = e_njk
    side(si) % b  = nv
  end if

  if(side(sj) % ea .eq. e_ijn) then
    side(sj) % ea = e_ink
    side(sj) % a  = nv
  end if
  if(side(sj) % eb .eq. e_ijn) then
    side(sj) % eb = e_ink
    side(sj) % b  = nv
  end if

  if(side(sk) % ea .eq. e_ijn) then
    side(sk) % ea = e_ijn
    side(sk) % a  = nv
  end if
  if(side(sk) % eb .eq. e_ijn) then
    side(sk) % eb = e_ijn
    side(sk) % b  = nv
  end if

! if(side(sj) % ea .eq. e_ijn) side(sj) % a = nv
! if(side(sj) % eb .eq. e_ijn) side(sj) % b = nv
! if(side(sk) % ea .eq. e_ijn) side(sk) % a = nv
! if(side(sk) % eb .eq. e_ijn) side(sk) % b = nv

!  ! One could think of this, but it's probably not needed
!  vert(i)  % nne = vert(i) % nne + 1
!  vert(j)  % nne = vert(j) % nne + 1
!  vert(k)  % nne = vert(k) % nne + 1
!  vert(nv) % nne = 3

  !--------------!
  !   Checking   !
  !--------------!
  do run = 1, 3
    if(run .eq. 1) s = si
    if(run .eq. 2) s = sj
    if(run .eq. 3) s = sk
    if(run .eq. 4) s = s_in
    if(run .eq. 5) s = s_jn
    if(run .eq. 6) s = s_kn
    ea = side(s) % ea
    if(ea > 0) then
      if( (elem(ea) % i + elem(ea) % j + elem(ea) % k) .ne.  &
          (side(s) % c + side(s) % d + side(s) % a) ) then
        print *, '# ERROR A in splitting an element at run', run
        stop
      end if
    end if

    eb = side(s) % eb
    if(eb > 0) then
      if( (elem(eb) % i + elem(eb) % j + elem(eb) % k) .ne.  &
          (side(s) % c + side(s) % d + side(s) % b) ) then
        print *, '# ERROR B in splitting an element at run', run
        stop
      end if
    end if
  end do

  do run = 1, 6
    if(run .eq. 1) e = e_ijn
    if(run .eq. 2) e = e_njk
    if(run .eq. 3) e = e_ink
    if(run .eq. 4) e = ei
    if(run .eq. 5) e = ej
    if(run .eq. 6) e = ek
    if(e > 0) then  ! ei, ej or ek can be zero
      sum_ijk = elem(e) % i + elem(e) % j + elem(e) % k
      sum_cd  = side(elem(e) % si) % c + side(elem(e) % si) % d  &
              + side(elem(e) % sj) % c + side(elem(e) % sj) % d  &
              + side(elem(e) % sk) % c + side(elem(e) % sk) % d
      if( sum_cd / sum_ijk .ne. 2 ) then
        print *, '# ERROR C in splitting an element!'
        stop
      end if
    end if
  end do

  end subroutine
