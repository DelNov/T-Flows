!==============================================================================!
  subroutine Front_Mod_Relax_Topology(front)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Front_Type), target :: front
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: vert(:)
  type(Side_Type), pointer :: side(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ns, ne
  integer                  :: a, b, c, d, s, t, e
!==============================================================================!

  ! Take aliases
  nv   => front % n_verts
  ns   => front % n_sides
  ne   => front % n_elems
  vert => front % vert
  side => front % side
  elem => front % elem

  call Front_Mod_Count_Vertex_Elements(front)
  call Front_Mod_Find_Boundaries(front)

  do t = 6, 3, -1
    do s = 1, ns

      a = side(s) % a
      b = side(s) % b
      c = side(s) % c
      d = side(s) % d

      ! This is how I check if side is on a boundary
      if( min(a, b, c, d) > 0 ) then

        if( .not. side(s) % boundary ) then

          e = vert(c) % nne + vert(d) % nne  &
            - vert(a) % nne - vert(b) % nne

          if(e .eq. t) then
            vert(a) % nne = vert(a) % nne + 1
            vert(b) % nne = vert(b) % nne + 1
            vert(c) % nne = vert(c) % nne - 1
            vert(d) % nne = vert(d) % nne - 1
            call Front_Mod_Swap_Side(front, s)
          end if

        end if
      end if  ! side is on a boundary

    end do
  end do

  end subroutine
