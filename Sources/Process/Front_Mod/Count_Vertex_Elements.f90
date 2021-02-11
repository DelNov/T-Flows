!==============================================================================!
  subroutine Front_Mod_Count_Vertex_Elements(front)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Front_Type), target :: front
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: vert(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ne
  integer                  :: e, i, j, k
!==============================================================================!

  ! Take aliases
  nv   => front % n_verts
  ne   => front % n_elems
  vert => front % vert
  elem => front % elem

  vert(1:nv) % nne = 0
  do e = 1, ne
    i = elem(e) % i
    j = elem(e) % j
    k = elem(e) % k
    vert(i) % nne = vert(i) % nne + 1
    vert(j) % nne = vert(j) % nne + 1
    vert(k) % nne = vert(k) % nne + 1
  end do

  end subroutine
