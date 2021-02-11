!==============================================================================!
  subroutine Surf_Mod_Count_Vertex_Elements(surf)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Surf_Type), target :: surf
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: vert(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ne
  integer                  :: e, i, j, k
!==============================================================================!

  ! Take aliases
  nv   => surf % n_verts
  ne   => surf % n_elems
  vert => surf % vert
  elem => surf % elem

  vert(1:nv) % nne = 0
  do e = 1, ne
    i = elem(e) % v(1)
    j = elem(e) % v(2)
    k = elem(e) % v(3)
    vert(i) % nne = vert(i) % nne + 1
    vert(j) % nne = vert(j) % nne + 1
    vert(k) % nne = vert(k) % nne + 1
  end do

  end subroutine
