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

  if(elem(ea) % ei .eq. eb) then
    ead = elem(ea) % ej;  eac = elem(ea) % ek
    sad = elem(ea) % sj;  sac = elem(ea) % sk
  end if
  if(elem(ea) % ej .eq. eb) then
    ead = elem(ea) % ek;  eac = elem(ea) % ei
    sad = elem(ea) % sk;  sac = elem(ea) % si
  end if
  if(elem(ea) % ek .eq. eb) then
    ead = elem(ea) % ei;  eac = elem(ea) % ej
    sad = elem(ea) % si;  sac = elem(ea) % sj
  end if

  if(elem(eb) % ei .eq. ea) then
    ebc = elem(eb) % ej;  ebd = elem(eb) % ek
    sbc = elem(eb) % sj;  sbd = elem(eb) % sk
  end if
  if(elem(eb) % ej .eq. ea) then
    ebc = elem(eb) % ek;  ebd = elem(eb) % ei
    sbc = elem(eb) % sk;  sbd = elem(eb) % si
  end if
  if(elem(eb) % ek .eq. ea) then
    ebc = elem(eb) % ei;  ebd = elem(eb) % ej
    sbc = elem(eb) % si;  sbd = elem(eb) % sj
  end if

  !--------------------------------------------!
  !   Change the orientation of the elements   !
  !--------------------------------------------!
  elem(ea) % i  = a;    elem(ea) % j  = b;    elem(ea) % k  = d
  elem(ea) % ei = ebd;  elem(ea) % ej = ead;  elem(ea) % ek = eb
  elem(ea) % si = sbd;  elem(ea) % sj = sad;  elem(ea) % sk = s

  elem(eb) % i  = a;    elem(eb) % j  = c;    elem(eb) % k  = b
  elem(eb) % ei = ebc;  elem(eb) % ej = ea;   elem(eb) % ek = eac
  elem(eb) % si = sbc;  elem(eb) % sj = s;    elem(eb) % sk = sac

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
    if(elem(eac) % ei .eq. ea) elem(eac) % ei = eb
    if(elem(eac) % ej .eq. ea) elem(eac) % ej = eb
    if(elem(eac) % ek .eq. ea) elem(eac) % ek = eb
  end if

  if(ebd .ne. 0) then
    if(elem(ebd) % ei .eq. eb) elem(ebd) % ei = ea
    if(elem(ebd) % ej .eq. eb) elem(ebd) % ej = ea
    if(elem(ebd) % ek .eq. eb) elem(ebd) % ek = ea
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
