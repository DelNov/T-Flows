!==============================================================================!
  subroutine Allocate_Surf(Surf, Flow)
!------------------------------------------------------------------------------!
!>  This subroutine is responsible for memory allocation and initialization of
!>  the surface object in Process. It assigns pointers to the relevant flow
!>  and grid objects.  The subroutine allocates memory for elements, vertices,
!>  and sides on the surface, along with data structures for data exchange if
!>  parallel runs.  The subroutine concludes by initializing the surface's
!>  local variables to set up for further processing.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf  !! parent class
  type(Field_Type), target :: Flow  !! flow field for which it is defined
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: Grid
!==============================================================================!

  ! Take aliases to object vertex Flow around
  Surf % pnt_flow => Flow
  Surf % pnt_grid => Flow % pnt_grid

  ! Take alias
  Grid => Flow % pnt_grid

  ! Surf shares surface mesh among processors
  Surf % mesh_divided = .false.

  ! Allocate memory
  allocate(Surf % Elem(Grid % n_cells))
  allocate(Surf % Vert(Grid % n_cells))
  allocate(Surf % side(Grid % n_cells * 3))

  ! Allocate logical array if cell holds vertices 
  ! (not sure if this will be needed)
  allocate(Surf % cell_has_vertex(Surf % pnt_grid % n_cells))
  Surf % cell_has_vertex(:) = .false.

  ! Allocate memory for buffers
  allocate(Surf % buff_x(Grid % n_cells));  Surf % buff_x(:) = 0.0
  allocate(Surf % buff_y(Grid % n_cells));  Surf % buff_y(:) = 0.0
  allocate(Surf % buff_z(Grid % n_cells));  Surf % buff_z(:) = 0.0
  allocate(Surf % buff_v(Grid % n_cells));  Surf % buff_v(:) = 0.0
  allocate(Surf % buff_i(Grid % n_cells));  Surf % buff_i(:) = 0
  allocate(Surf % buff_j(Grid % n_cells));  Surf % buff_j(:) = 0
  allocate(Surf % buff_k(Grid % n_cells));  Surf % buff_k(:) = 0
  allocate(Surf % buff_n(Grid % n_cells));  Surf % buff_n(:) = 0

  ! Initialize surf's local variables
  call Surf % Initialize_Surf()

  end subroutine
