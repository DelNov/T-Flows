!==============================================================================!
  subroutine Surf_Mod_Find_Vertex_Elements(surf)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Surf_Type), target :: surf
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: vert(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ne
  integer                  :: e, v, i, j, k
!==============================================================================!

  ! Take aliases
  nv   => surf % n_verts
  ne   => surf % n_elems
  vert => surf % vert
  elem => surf % elem

  call Surf_Mod_Count_Vertex_Elements(surf)

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
    i = elem(e) % v(1)
    j = elem(e) % v(2)
    k = elem(e) % v(3)
    vert(i) % nne = vert(i) % nne + 1;  vert(i) % vert_e(vert(i) % nne) = e
    vert(j) % nne = vert(j) % nne + 1;  vert(j) % vert_e(vert(j) % nne) = e
    vert(k) % nne = vert(k) % nne + 1;  vert(k) % vert_e(vert(k) % nne) = e
  end do

  end subroutine
