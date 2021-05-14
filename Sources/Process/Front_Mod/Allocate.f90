!==============================================================================!
  subroutine Allocate_Front(Front, flow)
!------------------------------------------------------------------------------!
!   Surface genesis                                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type), target :: Front
  type(Field_Type),  target :: flow
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: nb, nc, nf
!==============================================================================!

  ! Take aliases to object vertex flow around
  Front % pnt_flow => flow
  Front % pnt_grid => flow % pnt_grid

  ! Take aliases
  grid => flow % pnt_grid
  nb   =  grid % n_bnd_cells
  nc   =  grid % n_cells
  nf   =  grid % n_faces

  ! Allocate memory
  allocate(Front % elem(MAX_SURFACE_ELEMENTS))
  allocate(Front % vert(MAX_SURFACE_VERTICES))
  allocate(Front % side(MAX_SURFACE_ELEMENTS * 3))

  if(flow % mass_transfer) then
    allocate(Front % cell_at_elem(-nb:nc)); Front % cell_at_elem(-nb:nc) = 0
    allocate(Front % face_at_elem(  2,nf)); Front % face_at_elem(  :, :) = 0
  end if

  ! Initialize front's local variables
  call Front % Initialize_Front()

  end subroutine
