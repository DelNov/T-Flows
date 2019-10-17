!==============================================================================!
  subroutine Surf_Mod_Find_Boundaries(surf)
!------------------------------------------------------------------------------!
!   Looks for boundary sides and vertices                                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Surf_Type), target :: surf
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: vert(:)
  type(Side_Type), pointer :: side(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ns, ne
  integer                  :: s
!==============================================================================!

  ! Take aliases
  nv   => surf % n_verts
  ns   => surf % n_sides
  ne   => surf % n_elems
  vert => surf % vert
  side => surf % side
  elem => surf % elem

  ! First find boundary sides
  side(1:ns) % boundary = .false.
  do s = 1, ns
    if(min(side(s) % ea, side(s) % eb) .eq. 0) side(s) % boundary = .true.
  end do

  ! Then spread this information to nodes
  vert(1:nv) % boundary = .false.
  do s = 1, ns
    if(side(s) % boundary) then
      vert(side(s) % c) % boundary = .true.
      vert(side(s) % d) % boundary = .true.
    end if
  end do

  end subroutine
