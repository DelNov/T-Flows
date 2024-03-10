!==============================================================================!
  subroutine Find_Boundaries(Surf)
!------------------------------------------------------------------------------!
!>  Find_Boundaries identifies which sides and vertices of the surface mesh
!>  are located on the boundaries. This process involves checking each side
!>  to determine if it constitutes a boundary side and then propagating this
!>  information to the vertices connected to these sides.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf  !! parent class
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: Vert(:)
  type(Side_Type), pointer :: side(:)
  type(Elem_Type), pointer :: Elem(:)
  integer,         pointer :: nv, ns, ne
  integer                  :: s
!==============================================================================!

  ! Take aliases
  nv   => Surf % n_verts
  ns   => Surf % n_sides
  ne   => Surf % n_elems
  Vert => Surf % Vert
  side => Surf % side
  Elem => Surf % Elem

  ! First find boundary sides
  side(1:ns) % boundary = .false.
  do s = 1, ns
    if(min(side(s) % ea, side(s) % eb) .eq. 0) side(s) % boundary = .true.
  end do

  ! Then spread this information to vertices
  Vert(1:nv) % boundary = .false.
  do s = 1, ns
    if(side(s) % boundary) then
      Vert(side(s) % c) % boundary = .true.
      Vert(side(s) % d) % boundary = .true.
    end if
  end do

  end subroutine
