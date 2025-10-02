!==============================================================================!
  subroutine Relax_Geometry(Surf)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to optimize the geometric layout of a surface
!>  mesh. Unlike Relax_Topology, which uses the number of elements surrounding
!>  each node as a criterion for quality, Relax_Geometry evaluates and adjusts
!>  the mesh based on the ratio of the sides connecting nodes. It only leads
!>  to minor improvements in mesh qulity.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf  !! parent class
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: Vert(:)
  type(Side_Type), pointer :: side(:)
  type(Elem_Type), pointer :: Elem(:)
  integer,         pointer :: nv, ns, ne
  integer                  :: a, b, c, d, s, t
  real                     :: dist_ab, dist_cd
!==============================================================================!

  ! Take aliases
  nv   => Surf % n_verts
  ns   => Surf % n_sides
  ne   => Surf % n_elems
  Vert => Surf % Vert
  side => Surf % side
  Elem => Surf % Elem

  ! Boundary identification
  call Surf % Find_Boundaries()

  !-----------------------------------!
  !   Iterative geometry relaxation   !
  !- - - - - - - - - - - - - - - - - -+------------------------------!
  !   Executes a loop to adjust the geometry of the mesh, focusing   !
  !   on achieving a more uniform distribution of side lengths.      !
  !------------------------------------------------------------------!
  do t = 6, 3, -1
    do s = 1, ns

      a = side(s) % a
      b = side(s) % b
      c = side(s) % c
      d = side(s) % d

      ! This is how I check if side is on a boundary
      if( min(a, b, c, d) > 0 ) then

        if( .not. side(s) % boundary ) then

          dist_ab = Math % Distance(                         &
               Vert(a) % x_n, Vert(a) % y_n, Vert(a) % z_n,  &
               Vert(b) % x_n, Vert(b) % y_n, Vert(b) % z_n)
          dist_cd = Math % Distance(                         &
               Vert(c) % x_n, Vert(c) % y_n, Vert(c) % z_n,  &
               Vert(d) % x_n, Vert(d) % y_n, Vert(d) % z_n)

          if(dist_ab < dist_cd*TWO_THIRDS) then
            call Surf % Swap_Side(s)
          end if

        end if
      end if  ! side is on a boundary

    end do
  end do

  end subroutine
