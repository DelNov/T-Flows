!==============================================================================!
  subroutine Relax_Topology(Surf)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: vert(:)
  type(Side_Type), pointer :: side(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ns, ne
  integer                  :: a, b, c, d, s, t, e
!==============================================================================!

  ! Take aliases
  nv   => Surf % n_verts
  ns   => Surf % n_sides
  ne   => Surf % n_elems
  vert => Surf % vert
  side => Surf % side
  elem => Surf % elem

  call Surf % Find_Vertex_Elements()
  call Surf % Find_Boundaries()

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
            call Surf % Swap_Side(s)
          end if

        end if
      end if  ! side is on a boundary

    end do
  end do

  end subroutine
