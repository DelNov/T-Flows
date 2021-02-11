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
  integer                  :: e, v, i_v
!==============================================================================!

  ! Take aliases
  nv   => front % n_verts
  ne   => front % n_elems
  vert => front % vert
  elem => front % elem

  call Front_Mod_Count_Vertex_Elements(front)

  ! (Re)allocate memory for elements around each vertex
  do v = 1, nv
    if(allocated(vert(v) % vert_e)) then
      if(size(vert(v) % vert_e) .ne. vert(v) % nne) then
        deallocate(vert(v) % vert_e)
        allocate(vert(v) % vert_e(vert(v) % nne))
      end if
    else
      allocate(vert(v) % vert_e(vert(v) % nne))
    end if
  end do

  ! Store elements around each vertex
  vert(1:nv) % nne = 0
  do e = 1, ne
    do i_v = 1, elem(e) % nv
      v = elem(e) % v(i_v)
      vert(v) % nne = vert(v) % nne + 1;  vert(v) % vert_e(vert(v) % nne) = e
    end do
  end do

  end subroutine
