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

  ! Allocate logical array if cell holds vertices 
  ! (not sure if this will be needed)
  allocate(Surf % cell_has_vertex(Surf % pnt_grid % n_cells))
  Surf % cell_has_vertex(:) = .false.

  ! Allocate memory for working arrays
  ! (Not used yet, will be used in parallel version)
  ! allocate(i_work(Surf % n_verts * Surf % N_I_VARS))
  ! allocate(l_work(Surf % n_verts * Surf % N_L_VARS))
  ! allocate(r_work(Surf % n_verts * Surf % N_R_VARS))

  !-----------------------------!
  !   Initialize all elements   !
  !-----------------------------!
  allocate(Surf % elem(MAX_SURFACE_ELEMENTS))
  do e = 1, MAX_SURFACE_ELEMENTS
    Surf % elem(e) % nne  = 0
    Surf % elem(e) % v(:) = 0
    Surf % elem(e) % e(:) = 0
    Surf % elem(e) % s(:) = 0
  end do
  Surf % n_elems = 0

  !-----------------------------!
  !   Initialize all vertices   !
  !-----------------------------!
  allocate(Surf % vert(MAX_SURFACE_VERTICES))
  do v = 1, MAX_SURFACE_VERTICES

    Surf % vert(v) % nne = 0

    ! Set initial velocity to zero
    Surf % vert(v) % u_n = 0.0
    Surf % vert(v) % v_n = 0.0
    Surf % vert(v) % w_n = 0.0

    ! Set initial coordinates to zero
    Surf % vert(v) % x_n = 0.0
    Surf % vert(v) % y_n = 0.0
    Surf % vert(v) % z_n = 0.0

    Surf % vert(v) % x_o = 0.0
    Surf % vert(v) % y_o = 0.0
    Surf % vert(v) % z_o = 0.0

    ! Set initial cell, node and boundary cell to zero
    Surf % vert(v) % cell     = 0
    Surf % vert(v) % node     = 0
    Surf % vert(v) % bnd_cell = 0
    Surf % vert(v) % bnd_face = 0

    ! Assume vertex is in the domain
    ! (A smarter way could be worked out, depending ...
    ! ... on the result of the call to Find_Nearest_Cell)
    Surf % vert(v) % escaped   = .false.

    ! Is vertex in this processor?
    Surf % vert(v) % proc = 0
    Surf % vert(v) % buff = 0

  end do
  Surf % n_verts = 0

  !--------------------------!
  !   Initialize all sides   !
  !--------------------------!
  allocate(Surf % side(MAX_SURFACE_ELEMENTS * 3))
  do s = 1, MAX_SURFACE_ELEMENTS * 3
    Surf % side(s) % a        = 0
    Surf % side(s) % b        = 0
    Surf % side(s) % c        = 0
    Surf % side(s) % d        = 0
    Surf % side(s) % ei       = 0
    Surf % side(s) % ea       = 0
    Surf % side(s) % eb       = 0
    Surf % side(s) % length   = 0.0
    Surf % side(s) % boundary = .false.
  end do
  Surf % n_sides = 0

  end subroutine
