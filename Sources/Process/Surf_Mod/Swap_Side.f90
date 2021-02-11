!==============================================================================!
  subroutine Surf_Mod_Swap_Side(surf, s)
!------------------------------------------------------------------------------!
!  This function calculates radii of inscribed and circumscribed circle        !
!  for a given element (int e)                                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Surf_Type), target :: surf
  integer                 :: s
!-----------------------------------[Locals]-----------------------------------!
  integer                  :: a, b, c, d, ea, eb
  integer                  :: eac, ead, ebc, ebd
  integer                  :: sad, sac, sbc, sbd
  type(Vert_Type), pointer :: vert(:)
  type(Elem_Type), pointer :: elem(:)
  type(Side_Type), pointer :: side(:)
!==============================================================================!

  ! Take aliases
  vert => surf % vert
  elem => surf % elem
  side => surf % side

  ea = side(s) % ea
  eb = side(s) % eb
  a  = side(s) % a
  b  = side(s) % b
  c  = side(s) % c
  d  = side(s) % d

  eac = 0; ead = 0; ebc = 0; ebd = 0
  sad = 0; sac = 0; sbc = 0; sbd = 0

  if(elem(ea) % e(1) .eq. eb) then
    ead = elem(ea) % e(2);  eac = elem(ea) % e(3)
    sad = elem(ea) % s(2);  sac = elem(ea) % s(3)
  end if
  if(elem(ea) % e(2) .eq. eb) then
    ead = elem(ea) % e(3);  eac = elem(ea) % e(1)
    sad = elem(ea) % s(3);  sac = elem(ea) % s(1)
  end if
  if(elem(ea) % e(3) .eq. eb) then
    ead = elem(ea) % e(1);  eac = elem(ea) % e(2)
    sad = elem(ea) % s(1);  sac = elem(ea) % s(2)
  end if

  if(elem(eb) % e(1) .eq. ea) then
    ebc = elem(eb) % e(2);  ebd = elem(eb) % e(3)
    sbc = elem(eb) % s(2);  sbd = elem(eb) % s(3)
  end if
  if(elem(eb) % e(2) .eq. ea) then
    ebc = elem(eb) % e(3);  ebd = elem(eb) % e(1)
    sbc = elem(eb) % s(3);  sbd = elem(eb) % s(1)
  end if
  if(elem(eb) % e(3) .eq. ea) then
    ebc = elem(eb) % e(1);  ebd = elem(eb) % e(2)
    sbc = elem(eb) % s(1);  sbd = elem(eb) % s(2)
  end if

  !--------------------------------------------!
  !   Change the orientation of the elements   !
  !--------------------------------------------!
  elem(ea) % v(1) = a;    elem(ea) % v(2) = b;    elem(ea) % v(3) = d
  elem(ea) % e(1) = ebd;  elem(ea) % e(2) = ead;  elem(ea) % e(3) = eb
  elem(ea) % s(1) = sbd;  elem(ea) % s(2) = sad;  elem(ea) % s(3) = s

  elem(eb) % v(1) = a;    elem(eb) % v(2) = c;    elem(eb) % v(3) = b
  elem(eb) % e(1) = ebc;  elem(eb) % e(2) = ea;   elem(eb) % e(3) = eac
  elem(eb) % s(1) = sbc;  elem(eb) % s(2) = s;    elem(eb) % s(3) = sac

  !-------------------------------------------------!
  !   Change the orientation of the side properly   !
  !-------------------------------------------------!
  if(a < b) then
    side(s) % c  = a;   side(s) % d  = b
    side(s) % a  = d;   side(s) % b  = c
    side(s) % ea = ea;  side(s) % eb = eb
  else
    side(s) % c  = b;   side(s) % d  = a
    side(s) % a  = c;   side(s) % b  = d
    side(s) % ea = eb;  side(s) % eb = ea
  end if

  !-----------------------------------------------------!
  !   Update information on the neighbouring elements   !
  !-----------------------------------------------------!
  if(eac .ne. 0) then
    if(elem(eac) % e(1) .eq. ea) elem(eac) % e(1) = eb
    if(elem(eac) % e(2) .eq. ea) elem(eac) % e(2) = eb
    if(elem(eac) % e(3) .eq. ea) elem(eac) % e(3) = eb
  end if

  if(ebd .ne. 0) then
    if(elem(ebd) % e(1) .eq. eb) elem(ebd) % e(1) = ea
    if(elem(ebd) % e(2) .eq. eb) elem(ebd) % e(2) = ea
    if(elem(ebd) % e(3) .eq. eb) elem(ebd) % e(3) = ea
  end if

  !--------------------------------------------------!
  !   Update information on the neighbouring sides   !
  !--------------------------------------------------!
  if(side(sad) % ea .eq. ea) side(sad) % a = b
  if(side(sad) % eb .eq. ea) side(sad) % b = b

  if(side(sbc) % ea .eq. eb) side(sbc) % a = a
  if(side(sbc) % eb .eq. eb) side(sbc) % b = a

  if(side(sbd) % ea .eq. eb) then
    side(sbd) % ea = ea
    side(sbd) % a  = a
  end if
  if(side(sbd) % eb .eq. eb) then
    side(sbd) % eb = ea
    side(sbd) % b  = a
  end if

  if(side(sac) % ea .eq. ea) then
    side(sac) % ea = eb
    side(sac) % a  = b
  end if
  if(side(sac) % eb .eq. ea) then
    side(sac) % eb = eb
    side(sac) % b  = b
  end if

  end subroutine
