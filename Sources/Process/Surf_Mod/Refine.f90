!==============================================================================!
  subroutine Refine(Surf, n_biggest)
!------------------------------------------------------------------------------!
!   Refines ten biggest elements on the surface Surf                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf
  integer                  :: n_biggest  ! number of biggest elements for ref.
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: Vert(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ne
  real,    allocatable     :: areas(:)
  integer, allocatable     :: elems(:)
  integer                  :: ne_old, e
!==============================================================================!

  ! Take aliases
  nv   => Surf % n_verts
  ne   => Surf % n_elems
  Vert => Surf % Vert
  elem => Surf % elem

  allocate(areas(ne));  areas(:) = 0.0
  allocate(elems(ne));  elems(:) = 0

  do e = 1, ne
    elems(e) = e
    areas(e) = elem(e) % area
  end do

  ! Sort all elements by their areas
  call Sort % Real_Carry_Int(areas(1:ne), elems(1:ne))

  ! Refine ten biggest element
  ne_old = ne
  do e = ne_old, ne_old - n_biggest, -1
    call Surf % Split_Element(elems(e))
  end do

  end subroutine
