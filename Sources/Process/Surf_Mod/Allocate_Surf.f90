!==============================================================================!
  subroutine Allocate_Surf(Surf, Flow)
!------------------------------------------------------------------------------!
!   Surface genesis                                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf
  type(Field_Type), target :: Flow
!==============================================================================!

  ! Take aliases to object vertex Flow around
  Surf % pnt_flow => Flow
  Surf % pnt_grid => Flow % pnt_grid

  ! Surf shares surface mesh among processors
  Surf % mesh_divided = .false.

  ! Allocate memory
  allocate(Surf % Elem(MAX_SURFACE_ELEMENTS))
  allocate(Surf % Vert(MAX_SURFACE_VERTICES))
  allocate(Surf % side(MAX_SURFACE_ELEMENTS * 3))

  ! Allocate logical array if cell holds vertices 
  ! (not sure if this will be needed)
  allocate(Surf % cell_has_vertex(Surf % pnt_grid % n_cells))
  Surf % cell_has_vertex(:) = .false.

  ! Allocate memory for buffers
  allocate(Surf % buff_x(MAX_SURFACE_VERTICES));  Surf % buff_x(:) = 0.0
  allocate(Surf % buff_y(MAX_SURFACE_VERTICES));  Surf % buff_y(:) = 0.0
  allocate(Surf % buff_z(MAX_SURFACE_VERTICES));  Surf % buff_z(:) = 0.0
  allocate(Surf % buff_v(MAX_SURFACE_VERTICES));  Surf % buff_v(:) = 0.0
  allocate(Surf % buff_i(MAX_SURFACE_ELEMENTS));  Surf % buff_i(:) = 0
  allocate(Surf % buff_j(MAX_SURFACE_ELEMENTS));  Surf % buff_j(:) = 0
  allocate(Surf % buff_k(MAX_SURFACE_ELEMENTS));  Surf % buff_k(:) = 0
  allocate(Surf % buff_n(MAX_SURFACE_ELEMENTS));  Surf % buff_n(:) = 0

  ! Initialize surf's local variables
  call Surf % Initialize_Surf()

  end subroutine
