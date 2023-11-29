!==============================================================================!
  subroutine Allocate_Front(Front, Flow)
!------------------------------------------------------------------------------!
!>  This subroutine is responsible for memory allocation and initialization of
!>  the front structure in Process. It assigns pointers to the relevant flow
!>  and grid objects and sets the front to divide the surface mesh among
!>  processors. The subroutine allocates memory for elements, vertices, and
!>  sides of the front, along with boundary nodes and data structures for cell
!>  and face intersections.  The subroutine concludes by initializing the
!>  front's local variables to set up for further processing.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type), target :: Front  !! parent class
  type(Field_Type),  target :: Flow   !! flow field for which it is defined
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
  allocate(Front % intersects_face(nf));    Front % intersects_face(:) = .false.

  ! Face-based intersection with surface (needed for phase change)
  allocate(Front % xs(nf));  Front % xs(:) = 0.0
  allocate(Front % ys(nf));  Front % ys(:) = 0.0
  allocate(Front % zs(nf));  Front % zs(:) = 0.0

  ! Initialize front's local variables
  call Front % Initialize_Front()

  end subroutine
