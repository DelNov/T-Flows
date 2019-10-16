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
  integer :: v, e
!==============================================================================!

  ! Take aliases to object vertex flow around
  surf % pnt_flow => flow
  surf % pnt_grid => flow % pnt_grid

  ! Allocate logical array if cell holds vertices 
  ! (not sure if this will be needed)
  allocate(surf % cell_has_vertex(surf % pnt_grid % n_cells))
  surf % cell_has_vertex(:) = .false.

  ! Allocate memory for working arrays
  allocate(i_work(surf % n_verts * N_I_VARS))
  allocate(l_work(surf % n_verts * N_L_VARS))
  allocate(r_work(surf % n_verts * N_R_VARS))

  !-----------------------------!
  !   Initialize all elements   !
  !-----------------------------!
  do e = 1, surf % n_elems
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

  !-----------------------------!
  !   Initialize all vertices   !
  !-----------------------------!
  do v = 1, surf % n_verts

    surf % vert(e) % nne = 0

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

  end subroutine
