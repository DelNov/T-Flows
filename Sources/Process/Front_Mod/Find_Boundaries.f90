!==============================================================================!
  subroutine Front_Mod_Find_Boundaries(front)
!------------------------------------------------------------------------------!
!   Looks for boundary sides and vertices                                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Front_Type), target :: front
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: vert(:)
  type(Side_Type), pointer :: side(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ns, ne
  integer                  :: s
!==============================================================================!

  ! Take aliases
  nv   => front % n_verts
  ns   => front % n_sides
  ne   => front % n_elems
  vert => front % vert
  side => front % side
  elem => front % elem

  ! First find boundary sides
  side(1:ns) % boundary = .false.
  do s = 1, ns
    if(min(side(s) % ea, side(s) % eb) .eq. 0) side(s) % boundary = .true.
  end do

  ! Then spread this information to vertices
  vert(1:nv) % boundary = .false.
  do s = 1, ns
    if(side(s) % boundary) then
      vert(side(s) % c) % boundary = .true.
      vert(side(s) % d) % boundary = .true.
    end if
  end do

  end subroutine
