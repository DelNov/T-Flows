!==============================================================================!
  subroutine Find_Vertex_Elements(Front)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type), target :: Front
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: vert(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ne
  integer                  :: e, v, i_ver
!==============================================================================!

  ! Take aliases
  nv   => Front % n_verts
  ne   => Front % n_elems
  vert => Front % vert
  elem => Front % elem

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
