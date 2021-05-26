!==============================================================================!
  subroutine Relax_Topology(Surf)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: Vert(:)
  type(Side_Type), pointer :: side(:)
  type(Elem_Type), pointer :: Elem(:)
  integer,         pointer :: nv, ns, ne
  integer                  :: a, b, c, d, s, t, e
!==============================================================================!

  ! Take aliases
  nv   => Surf % n_verts
  ns   => Surf % n_sides
  ne   => Surf % n_elems
  Vert => Surf % Vert
  side => Surf % side
  Elem => Surf % Elem

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

          e = Vert(c) % nne + Vert(d) % nne  &
            - Vert(a) % nne - Vert(b) % nne

          if(e .eq. t) then
            Vert(a) % nne = Vert(a) % nne + 1
            Vert(b) % nne = Vert(b) % nne + 1
            Vert(c) % nne = Vert(c) % nne - 1
            Vert(d) % nne = Vert(d) % nne - 1
            call Surf % Swap_Side(s)
          end if

        end if
      end if  ! side is on a boundary

    end do
  end do

  end subroutine
