!==============================================================================!
  subroutine Allocate_Front(Front, Flow)
!------------------------------------------------------------------------------!
!   Memory allocation for front                                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type), target :: Front
  type(Field_Type),  target :: Flow
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: nb, nc, nf
!==============================================================================!

  ! Take aliases to object vertex Flow around
  Front % pnt_flow => Flow
  Front % pnt_grid => Flow % pnt_grid

  ! Use just Grid from now on
  Grid => Flow % pnt_grid

  ! Front divides the surface mesh among processors
  Front % mesh_divided = .true.

  ! Allocate memory
  allocate(Front % Elem(max(Grid % n_cells,    &
                            Grid % n_nodes)))
  allocate(Front % Vert(max(Grid % n_cells,    &
                            Grid % n_nodes)))
  allocate(Front % side(max(Grid % n_cells,    &
                            Grid % n_nodes) * 3))

  allocate(Front % b_node_1(max(Grid % n_cells, Grid % n_nodes)));
  allocate(Front % b_node_2(max(Grid % n_cells, Grid % n_nodes)));
  Front % b_node_1(:) = 0
  Front % b_node_2(:) = 0

  nb =  Grid % n_bnd_cells
  nc =  Grid % n_cells
  nf =  Grid % n_faces
  allocate(Front % elem_in_cell (-nb:nc));  Front % elem_in_cell (:)   = 0
  allocate(Front % elems_at_face  (2,nf));  Front % elems_at_face(:,:) = 0

  ! Initialize front's local variables
  call Front % Initialize_Front()

  end subroutine
