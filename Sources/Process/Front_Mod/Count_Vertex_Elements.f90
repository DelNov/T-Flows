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
  integer                  :: e, v, i_v
!==============================================================================!

  ! Take aliases
  nv   => front % n_verts
  ne   => front % n_elems
  vert => front % vert
  elem => front % elem

  vert(1:nv) % nne = 0
  do e = 1, ne
    do i_v = 1, elem(e) % nv
      v = elem(e) % v(i_v)
      vert(v) % nne = vert(v) % nne + 1
    end do
  end do

  end subroutine
