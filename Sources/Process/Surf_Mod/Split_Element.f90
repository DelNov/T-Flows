!==============================================================================!
  subroutine Split_Element(Surf, e)
!------------------------------------------------------------------------------!
!   Create new vertex in the center of existing element specified with e       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf
  integer                  :: e
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: Vert(:)
  type(Side_Type), pointer :: side(:)
  type(Elem_Type), pointer :: Elem(:)
  integer,         pointer :: nv, ns, ne
  integer                  :: i, j, k, ei, ej, ek, s, si, sj, sk
  integer                  :: sum_ijk, sum_cd, run, ea, eb
  integer                  :: s_in, s_jn, s_kn, e_ijn, e_njk, e_ink
!==============================================================================!

  ! Take aliases
  nv   => Surf % n_verts
  ns   => Surf % n_sides
  ne   => Surf % n_elems
  Vert => Surf % Vert
  side => Surf % side
  Elem => Surf % Elem

  i = Elem(e) % v(1)
  j = Elem(e) % v(2)
  k = Elem(e) % v(3)

  ei = Elem(e) % e(1)
  ej = Elem(e) % e(2)
  ek = Elem(e) % e(3)

  si = Elem(e) % s(1)
  sj = Elem(e) % s(2)
  sk = Elem(e) % s(3)

  !-----------------------------------------------!
  !   Create new vertex in the element's center   !
  !-----------------------------------------------!
  nv = nv + 1
  Vert(nv) % x_n = (Vert(i) % x_n + Vert(j) % x_n + Vert(k) % x_n) * ONE_THIRD
  Vert(nv) % y_n = (Vert(i) % y_n + Vert(j) % y_n + Vert(k) % y_n) * ONE_THIRD
  Vert(nv) % z_n = (Vert(i) % z_n + Vert(j) % z_n + Vert(k) % z_n) * ONE_THIRD

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

  ! Handle new elements' vertices
  Elem(e_ijn) % v(1) = i
  Elem(e_ijn) % v(2) = j
  Elem(e_ijn) % v(3) = nv

  Elem(e_njk) % v(1) = nv
  Elem(e_njk) % v(2) = j
  Elem(e_njk) % v(3) = k

  Elem(e_ink) % v(1) = i
  Elem(e_ink) % v(2) = nv
  Elem(e_ink) % v(3) = k

  ! Handle new elements' neighbours
  Elem(e_ijn) % e(1) = e_njk
  Elem(e_ijn) % e(2) = e_ink
  Elem(e_ijn) % e(3) = ek

  Elem(e_njk) % e(1) = ei
  Elem(e_njk) % e(2) = e_ink
  Elem(e_njk) % e(3) = e_ijn

  Elem(e_ink) % e(1) = e_njk
  Elem(e_ink) % e(2) = ej
  Elem(e_ink) % e(3) = e_ijn

  ! Handle new elements' sides
  Elem(e_ijn) % s(1) = s_jn
  Elem(e_ijn) % s(2) = s_in
  Elem(e_ijn) % s(3) = sk

  Elem(e_njk) % s(1) = si
  Elem(e_njk) % s(2) = s_kn
  Elem(e_njk) % s(3) = s_jn

  Elem(e_ink) % s(1) = s_kn
  Elem(e_ink) % s(2) = sj
  Elem(e_ink) % s(3) = s_in

  !-------------------------------------------!
  !   Increase number of sides and elements   !
  !-------------------------------------------!
  ns = ns + 3
  ne = ne + 2

  !-------------------------------------!
  !   Correct neighbours where needed   !
  !-------------------------------------!

  if(ei .gt. 0) then
    if(Elem(ei) % e(1) .eq. e_ijn) Elem(ei) % e(1) = e_njk
    if(Elem(ei) % e(2) .eq. e_ijn) Elem(ei) % e(2) = e_njk
    if(Elem(ei) % e(3) .eq. e_ijn) Elem(ei) % e(3) = e_njk
  end if

  if(ej .gt. 0) then
    if(Elem(ej) % e(1) .eq. e_ijn) Elem(ej) % e(1) = e_ink
    if(Elem(ej) % e(2) .eq. e_ijn) Elem(ej) % e(2) = e_ink
    if(Elem(ej) % e(3) .eq. e_ijn) Elem(ej) % e(3) = e_ink
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
!  Vert(i)  % nne = Vert(i) % nne + 1
!  Vert(j)  % nne = Vert(j) % nne + 1
!  Vert(k)  % nne = Vert(k) % nne + 1
!  Vert(nv) % nne = 3

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
      if( (Elem(ea) % v(1) + Elem(ea) % v(2) + Elem(ea) % v(3)) .ne.  &
          (side(s) % c + side(s) % d + side(s) % a) ) then
        print *, '# ERROR A in splitting an element at run', run
        stop
      end if
    end if

    eb = side(s) % eb
    if(eb > 0) then
      if( (Elem(eb) % v(1) + Elem(eb) % v(2) + Elem(eb) % v(3)) .ne.  &
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
      sum_ijk = Elem(e) % v(1) + Elem(e) % v(2) + Elem(e) % v(3)
      sum_cd  = side(Elem(e) % s(1)) % c + side(Elem(e) % s(1)) % d  &
              + side(Elem(e) % s(2)) % c + side(Elem(e) % s(2)) % d  &
              + side(Elem(e) % s(3)) % c + side(Elem(e) % s(3)) % d
      if( sum_cd / sum_ijk .ne. 2 ) then
        print *, '# ERROR C in splitting an element!'
        stop
      end if
    end if
  end do

  end subroutine
