!==============================================================================!
  subroutine Swap_Side(Surf, s)
!------------------------------------------------------------------------------!
!>  Swap_Side focuses on reconfiguring the mesh by altering the connections
!>  and orientations of elements and sides around a specified side (s). This
!>  process is crucial for implementation of the algorithm in Relax_Topology
!>  and Relax_Geometry. This algorithm of mesh relaxation is taken from
!>  TRIPOS (https://github.com/Niceno/TRIPOS)
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Alias setup:                                                             !
!     - Establishes aliases for pointers to vertices (Vert), elements (Elem),  !
!      and sides (side) to simplify code navigation and enhance clarity.       !
!   * Retrieving side and element information:                                 !
!     - Retrieves information about the target side (s) and the elements (ea   !
!       and eb) adjacent to this side. Also, obtains the vertices (a, b, c, d) !
!       that are connected by this side and its adjacent elements.             !
!   * Initialization of variables:                                             !
!     - Initializes variables to track the neighboring elements and sides of   !
!       ea and eb.                                                             !
!   * Reconfiguring elements and sides:                                        !
!     - Changes the orientation of elements ea and eb by adjusting their       !
!       vertex and side connections. This reconfiguration is based on the      !
!       current mesh structure and the intended changes in the topology.       !
!     - Adjusts the orientation of the target side (s) to align with the       !
!       new configuration of elements ea and eb.                               !
!   * Updating neighboring elements and sides:                                 !
!     - Updates the neighboring elements and sides of ea and eb to reflect     !
!       the new mesh configuration. This step ensures consistency in the       !
!       mesh's topological structure after the swapping operation.             !
!   * Consistency checks:                                                      !
!     - Performs consistency checks to verify that the new topology is valid   !
!       and does not introduce any errors or inconsistencies in the mesh.      !
!       This is a critical step to maintain the integrity of the mesh.         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf  !! parent class
  integer                  :: s     !! side to swap
!-----------------------------------[Locals]-----------------------------------!
  integer                  :: a, b, c, d, ea, eb
  integer                  :: eac, ead, ebc, ebd
  integer                  :: sad, sac, sbc, sbd
  type(Vert_Type), pointer :: Vert(:)
  type(Elem_Type), pointer :: Elem(:)
  type(Side_Type), pointer :: side(:)
!==============================================================================!

  ! Take aliases
  Vert => Surf % Vert
  Elem => Surf % Elem
  side => Surf % side

  ea = side(s) % ea
  eb = side(s) % eb
  a  = side(s) % a
  b  = side(s) % b
  c  = side(s) % c
  d  = side(s) % d

  eac = 0; ead = 0; ebc = 0; ebd = 0
  sad = 0; sac = 0; sbc = 0; sbd = 0

  if(Elem(ea) % e(1) .eq. eb) then
    ead = Elem(ea) % e(2);  eac = Elem(ea) % e(3)
    sad = Elem(ea) % s(2);  sac = Elem(ea) % s(3)
  end if
  if(Elem(ea) % e(2) .eq. eb) then
    ead = Elem(ea) % e(3);  eac = Elem(ea) % e(1)
    sad = Elem(ea) % s(3);  sac = Elem(ea) % s(1)
  end if
  if(Elem(ea) % e(3) .eq. eb) then
    ead = Elem(ea) % e(1);  eac = Elem(ea) % e(2)
    sad = Elem(ea) % s(1);  sac = Elem(ea) % s(2)
  end if

  if(Elem(eb) % e(1) .eq. ea) then
    ebc = Elem(eb) % e(2);  ebd = Elem(eb) % e(3)
    sbc = Elem(eb) % s(2);  sbd = Elem(eb) % s(3)
  end if
  if(Elem(eb) % e(2) .eq. ea) then
    ebc = Elem(eb) % e(3);  ebd = Elem(eb) % e(1)
    sbc = Elem(eb) % s(3);  sbd = Elem(eb) % s(1)
  end if
  if(Elem(eb) % e(3) .eq. ea) then
    ebc = Elem(eb) % e(1);  ebd = Elem(eb) % e(2)
    sbc = Elem(eb) % s(1);  sbd = Elem(eb) % s(2)
  end if

  !--------------------------------------------!
  !   Change the orientation of the elements   !
  !--------------------------------------------!
  Elem(ea) % v(1) = a;    Elem(ea) % v(2) = b;    Elem(ea) % v(3) = d
  Elem(ea) % e(1) = ebd;  Elem(ea) % e(2) = ead;  Elem(ea) % e(3) = eb
  Elem(ea) % s(1) = sbd;  Elem(ea) % s(2) = sad;  Elem(ea) % s(3) = s

  Elem(eb) % v(1) = a;    Elem(eb) % v(2) = c;    Elem(eb) % v(3) = b
  Elem(eb) % e(1) = ebc;  Elem(eb) % e(2) = ea;   Elem(eb) % e(3) = eac
  Elem(eb) % s(1) = sbc;  Elem(eb) % s(2) = s;    Elem(eb) % s(3) = sac

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
    if(Elem(eac) % e(1) .eq. ea) Elem(eac) % e(1) = eb
    if(Elem(eac) % e(2) .eq. ea) Elem(eac) % e(2) = eb
    if(Elem(eac) % e(3) .eq. ea) Elem(eac) % e(3) = eb
  end if

  if(ebd .ne. 0) then
    if(Elem(ebd) % e(1) .eq. eb) Elem(ebd) % e(1) = ea
    if(Elem(ebd) % e(2) .eq. eb) Elem(ebd) % e(2) = ea
    if(Elem(ebd) % e(3) .eq. eb) Elem(ebd) % e(3) = ea
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
