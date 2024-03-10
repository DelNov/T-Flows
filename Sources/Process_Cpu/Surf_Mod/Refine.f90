!==============================================================================!
  subroutine Refine(Surf, n_biggest)
!------------------------------------------------------------------------------!
!>  This subroutine is tasked with refining the surface mesh by selectively
!>  targeting and subdividing the largest elements. By specifying the number
!>  of elements to refine (through n_biggest), the subroutine can concentrate
!>  refinement efforts on the most significant elements.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf       !! parent class
  integer                  :: n_biggest  !! number of biggest elements
                                         !! for refinement
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: Vert(:)
  type(Elem_Type), pointer :: Elem(:)
  integer,         pointer :: nv, ne
  real,    allocatable     :: areas(:)
  integer, allocatable     :: elems(:)
  integer                  :: ne_old, e
!==============================================================================!

  ! Take aliases
  nv   => Surf % n_verts
  ne   => Surf % n_elems
  Vert => Surf % Vert
  Elem => Surf % Elem

  allocate(areas(ne));  areas(:) = 0.0
  allocate(elems(ne));  elems(:) = 0

  do e = 1, ne
    elems(e) = e
    areas(e) = Elem(e) % area
  end do

  ! Sort all elements by their areas
  call Sort % Real_Carry_Int(areas(1:ne), elems(1:ne))

  ! Refine ten biggest element
  ne_old = ne
  do e = ne_old, ne_old - n_biggest, -1
    call Surf % Split_Element(elems(e))
  end do

  end subroutine
