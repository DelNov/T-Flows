!==============================================================================!
  subroutine Find_Vertex_Elements(Front)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type), target :: Front
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: Vert(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ne
  integer                  :: e, v, i_ver
!==============================================================================!

  ! Take aliases
  nv   => Front % n_verts
  ne   => Front % n_elems
  Vert => Front % Vert
  elem => Front % elem

  ! Initialize to zero
  do v = 1, nv
    Vert(v) % nne  = 0
    Vert(v) % e(:) = 0
  end do

  ! Store elements around each vertex (no checking!)
  do e = 1, ne
    do i_ver = 1, elem(e) % nv
      v = elem(e) % v(i_ver)
      Vert(v) % nne = Vert(v) % nne + 1;  Vert(v) % e(Vert(v) % nne) = e
    end do
  end do

  end subroutine
