!==============================================================================!
  subroutine Relax_Geometry(Surf)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: vert(:)
  type(Side_Type), pointer :: side(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ns, ne
  integer                  :: a, b, c, d, s, t
  real                     :: dist_ab, dist_cd
!==============================================================================!

  ! Take aliases
  nv   => Surf % n_verts
  ns   => Surf % n_sides
  ne   => Surf % n_elems
  vert => Surf % vert
  side => Surf % side
  elem => Surf % elem

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

          dist_ab = Math_Mod_Distance(                       &
               vert(a) % x_n, vert(a) % y_n, vert(a) % z_n,  &
               vert(b) % x_n, vert(b) % y_n, vert(b) % z_n)
          dist_cd = Math_Mod_Distance(                       &
               vert(c) % x_n, vert(c) % y_n, vert(c) % z_n,  &
               vert(d) % x_n, vert(d) % y_n, vert(d) % z_n)

          if(dist_ab < dist_cd*TWO_THIRDS) then
            call Surf % Swap_Side(s)
          end if

        end if
      end if  ! side is on a boundary

    end do
  end do

  end subroutine
