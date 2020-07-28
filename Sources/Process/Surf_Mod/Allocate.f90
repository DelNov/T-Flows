!==============================================================================!
  subroutine Surf_Mod_Allocate(surf, flow)
!------------------------------------------------------------------------------!
!   Surface genesis                                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Surf_Type),  target :: surf
  type(Field_Type), target :: flow
!----------------------------------[Locals]------------------------------------!
  integer :: v, e, s
!==============================================================================!

  ! Take aliases to object vertex flow around
  surf % pnt_flow => flow
  surf % pnt_grid => flow % pnt_grid

  ! Allocate logical array if cell holds vertices 
  ! (not sure if this will be needed)
  allocate(surf % cell_has_vertex(surf % pnt_grid % n_cells))
  surf % cell_has_vertex(:) = .false.

  ! Allocate memory for working arrays
  ! (Not used yet, will be used in parallel version)
  ! allocate(i_work(surf % n_verts * surf % N_I_VARS))
  ! allocate(l_work(surf % n_verts * surf % N_L_VARS))
  ! allocate(r_work(surf % n_verts * surf % N_R_VARS))

  !-----------------------------!
  !   Initialize all elements   !
  !-----------------------------!
  allocate(surf % elem(MAX_SURFACE_ELEMENTS))
  do e = 1, MAX_SURFACE_ELEMENTS
    surf % elem(e) % nne = 0
    surf % elem(e) %   i = 0
    surf % elem(e) %   j = 0
    surf % elem(e) %   k = 0
    surf % elem(e) %  ei = 0
    surf % elem(e) %  ej = 0
    surf % elem(e) %  ek = 0
    surf % elem(e) %  si = 0
    surf % elem(e) %  sj = 0
    surf % elem(e) %  sk = 0
  end do
  surf % n_elems = 0

  !-----------------------------!
  !   Initialize all vertices   !
  !-----------------------------!
  allocate(surf % vert(MAX_SURFACE_VERTICES))
  do v = 1, MAX_SURFACE_VERTICES

    surf % vert(v) % nne = 0

    ! Set initial velocity to zero
    surf % vert(v) % u = 0.0
    surf % vert(v) % v = 0.0
    surf % vert(v) % w = 0.0

    ! Set initial coordinates to zero
    surf % vert(v) % x_n = 0.0
    surf % vert(v) % y_n = 0.0
    surf % vert(v) % z_n = 0.0

    surf % vert(v) % x_o = 0.0
    surf % vert(v) % y_o = 0.0
    surf % vert(v) % z_o = 0.0

    ! Set initial cell, node and boundary cell to zero
    surf % vert(v) % cell     = 0
    surf % vert(v) % node     = 0
    surf % vert(v) % bnd_cell = 0
    surf % vert(v) % bnd_face = 0

    ! Assume vertex is in the domain
    ! (A smarter way could be worked out, depending ...
    ! ... on the result of the call to Find_Nearest_Cell)
    surf % vert(v) % escaped   = .false.

    ! Is vertex in this processor?
    surf % vert(v) % proc = 0
    surf % vert(v) % buff = 0

  end do
  surf % n_verts = 0

  !--------------------------!
  !   Initialize all sides   !
  !--------------------------!
  allocate(surf % side(MAX_SURFACE_ELEMENTS * 3))
  do s = 1, MAX_SURFACE_ELEMENTS * 3
    surf % side(s) % a        = 0
    surf % side(s) % b        = 0
    surf % side(s) % c        = 0
    surf % side(s) % d        = 0
    surf % side(s) % ei       = 0
    surf % side(s) % ea       = 0
    surf % side(s) % eb       = 0
    surf % side(s) % length   = 0.0
    surf % side(s) % boundary = .false.
  end do
  surf % n_sides = 0

  end subroutine
