!==============================================================================!
  subroutine Front_Mod_Allocate(front, flow)
!------------------------------------------------------------------------------!
!   Surface genesis                                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Front_Type), target :: front
  type(Field_Type), target :: flow
!----------------------------------[Locals]------------------------------------!
  integer :: v, e, s
!==============================================================================!

  ! Take aliases to object vertex flow around
  front % pnt_flow => flow
  front % pnt_grid => flow % pnt_grid

  ! Allocate logical array if cell holds vertices 
  ! (not sure if this will be needed)
  allocate(front % cell_has_vertex(front % pnt_grid % n_cells))
  front % cell_has_vertex(:) = .false.

  ! Allocate memory for working arrays
  ! (Not used yet, will be used in parallel version)
  ! allocate(i_work(front % n_verts * front % N_I_VARS))
  ! allocate(l_work(front % n_verts * front % N_L_VARS))
  ! allocate(r_work(front % n_verts * front % N_R_VARS))

  !-----------------------------!
  !   Initialize all elements   !
  !-----------------------------!
  allocate(front % elem(MAX_SURFACE_ELEMENTS))
  do e = 1, MAX_SURFACE_ELEMENTS
    front % elem(e) % nne = 0
    front % elem(e) %   i = 0
    front % elem(e) %   j = 0
    front % elem(e) %   k = 0
    front % elem(e) %  ei = 0
    front % elem(e) %  ej = 0
    front % elem(e) %  ek = 0
    front % elem(e) %  si = 0
    front % elem(e) %  sj = 0
    front % elem(e) %  sk = 0
  end do
  front % n_elems = 0

  !-----------------------------!
  !   Initialize all vertices   !
  !-----------------------------!
  allocate(front % vert(MAX_SURFACE_VERTICES))
  do v = 1, MAX_SURFACE_VERTICES

    front % vert(v) % nne = 0

    ! Set initial velocity to zero
    front % vert(v) % u = 0.0
    front % vert(v) % v = 0.0
    front % vert(v) % w = 0.0

    ! Set initial coordinates to zero
    front % vert(v) % x_n = 0.0
    front % vert(v) % y_n = 0.0
    front % vert(v) % z_n = 0.0

    front % vert(v) % x_o = 0.0
    front % vert(v) % y_o = 0.0
    front % vert(v) % z_o = 0.0

    ! Set initial cell, node and boundary cell to zero
    front % vert(v) % cell     = 0
    front % vert(v) % node     = 0
    front % vert(v) % bnd_cell = 0
    front % vert(v) % bnd_face = 0

    ! Assume vertex is in the domain
    ! (A smarter way could be worked out, depending ...
    ! ... on the result of the call to Find_Nearest_Cell)
    front % vert(v) % escaped   = .false.

    ! Is vertex in this processor?
    front % vert(v) % proc = 0
    front % vert(v) % buff = 0

  end do
  front % n_verts = 0

  !--------------------------!
  !   Initialize all sides   !
  !--------------------------!
  allocate(front % side(MAX_SURFACE_ELEMENTS * 3))
  do s = 1, MAX_SURFACE_ELEMENTS * 3
    front % side(s) % a        = 0
    front % side(s) % b        = 0
    front % side(s) % c        = 0
    front % side(s) % d        = 0
    front % side(s) % ei       = 0
    front % side(s) % ea       = 0
    front % side(s) % eb       = 0
    front % side(s) % length   = 0.0
    front % side(s) % boundary = .false.
  end do
  front % n_sides = 0

  end subroutine
