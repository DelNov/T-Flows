!==============================================================================!
  subroutine Front_Mod_Find_Vertex_Elements(front)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Front_Type), target :: front
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: vert(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ne
  integer                  :: e, v, i_ver
!==============================================================================!

  ! Take aliases
  nv   => front % n_verts
  ne   => front % n_elems
  vert => front % vert
  elem => front % elem

  ! Initialize to zero
  do v = 1, nv
    vert(v) % nne  = 0
    vert(v) % e(:) = 0
  end do

  ! Store elements around each vertex (no checking!)
  do e = 1, ne
    do i_ver = 1, elem(e) % nv
      v = elem(e) % v(i_ver)
      vert(v) % nne = vert(v) % nne + 1;  vert(v) % e(vert(v) % nne) = e
    end do
  end do

  end subroutine
