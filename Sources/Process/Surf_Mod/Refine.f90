!==============================================================================!
  subroutine Surf_Mod_Refine(surf, n_biggest)
!------------------------------------------------------------------------------!
!   Refines ten biggest elements on the surface surf                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Surf_Type), target :: surf
  integer                 :: n_biggest  ! number of biggest elements for ref.
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: vert(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ne
  real,    allocatable     :: areas(:)
  integer, allocatable     :: elems(:)
  integer                  :: ne_old, e
!==============================================================================!

  ! Take aliases
  nv   => surf % n_verts
  ne   => surf % n_elems
  vert => surf % vert
  elem => surf % elem

  allocate(areas(ne));  areas(:) = 0.0
  allocate(elems(ne));  elems(:) = 0

  do e = 1, ne
    elems(e) = e
    areas(e) = elem(e) % area
  end do

  ! Sort all elements by their areas
  call Sort_Mod_Real_Carry_Int(areas(1:ne), elems(1:ne))

  ! Refine ten biggest element
  ne_old = ne
  do e = ne_old, ne_old - n_biggest, -1
    call Surf_Mod_Split_Element(surf, elems(e))
  end do

  end subroutine
