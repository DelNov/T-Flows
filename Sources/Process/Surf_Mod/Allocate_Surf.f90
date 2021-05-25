!==============================================================================!
  subroutine Allocate_Surf(Surf, Flow)
!------------------------------------------------------------------------------!
!   Surface genesis                                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf
  type(Field_Type), target :: Flow
!----------------------------------[Locals]------------------------------------!
  integer :: v, e, s
!==============================================================================!

  ! Take aliases to object vertex Flow around
  Surf % pnt_flow => Flow
  Surf % pnt_grid => Flow % pnt_grid

  ! Allocate memory
  allocate(Surf % elem(MAX_SURFACE_ELEMENTS))
  allocate(Surf % Vert(MAX_SURFACE_VERTICES))
  allocate(Surf % side(MAX_SURFACE_ELEMENTS * 3))

  ! Allocate logical array if cell holds vertices 
  ! (not sure if this will be needed)
  allocate(Surf % cell_has_vertex(Surf % pnt_grid % n_cells))
  Surf % cell_has_vertex(:) = .false.

  ! Allocate memory for working arrays
  ! (Not used yet, will be used in parallel version)
  ! allocate(i_work(Surf % n_verts * Surf % N_I_VARS))
  ! allocate(l_work(Surf % n_verts * Surf % N_L_VARS))
  ! allocate(r_work(Surf % n_verts * Surf % N_R_VARS))

  ! Initialize surf's local variables
  call Surf % Initialize_Surf()

  end subroutine
