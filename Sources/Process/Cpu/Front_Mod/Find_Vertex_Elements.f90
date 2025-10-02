!==============================================================================!
  subroutine Find_Vertex_Elements(Front)
!------------------------------------------------------------------------------!
!>  This subroutine plays a crucial role in determining the relationship
!>  between vertices and elements of a front.  It is typically called before
!>  executing Calculate_Element_Centroids and Calculate_Element_Normals, which
!>  rely on the vertex-to-element relationship established here.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initializing each vertex's element count and associated elements to zero.!
!   * Iterating over each element to identify the vertices it contains.        !
!   * For each vertex, incrementing the number of neighboring elements and     !
!     storing the element identifiers.                                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type), target :: Front  !! parent class
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: Vert(:)
  type(Elem_Type), pointer :: Elem(:)
  integer,         pointer :: nv, ne
  integer                  :: e, v, i_ver
!==============================================================================!

  ! Take aliases
  nv   => Front % n_verts
  ne   => Front % n_elems
  Vert => Front % Vert
  Elem => Front % Elem

  ! Initialize to zero
  do v = 1, nv
    Vert(v) % nne  = 0
    Vert(v) % e(:) = 0
  end do

  ! Store elements around each vertex (no checking!)
  do e = 1, ne
    do i_ver = 1, Elem(e) % nv
      v = Elem(e) % v(i_ver)
      Vert(v) % nne = Vert(v) % nne + 1;  Vert(v) % e(Vert(v) % nne) = e
    end do
  end do

  end subroutine
