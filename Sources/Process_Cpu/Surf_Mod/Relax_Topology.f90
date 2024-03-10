!==============================================================================!
  subroutine Relax_Topology(Surf)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to modify the topological arrangement of the
!>  surface mesh elements and sides in T-Flows. It focuses on redistributing
!>  the mesh elements and adjusting their connections to maximize the number
!>  of nodes surrounded with six elements. This algorithm of mesh relaxation
!>  is taken from TRIPOS (https://github.com/Niceno/TRIPOS)
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf  !! parent class
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

  ! Update the mesh's current topological information, such as
  ! vertex-to-element connections and boundary identification.
  call Surf % Find_Vertex_Elements()
  call Surf % Find_Boundaries()

  !-----------------------------------!
  !   Iterative topology relaxation   !
  !- - - - - - - - - - - - - - - - - -+-------------------------------------+
  !   Executes a loop to relax the topology, starting with higher-valence   !
  !   vertices (vertices connected to more elements) and moving towards     !
  !   lower-valence ones. This gradual approach helps to uniformly          !
  !   distribute topological adjustments across the mesh.                   !
  !-------------------------------------------------------------------------!
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
