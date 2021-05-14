!==============================================================================!
  subroutine Initialize_Front(Front)
!------------------------------------------------------------------------------!
!   Surface genesis                                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type), target :: Front
!----------------------------------[Locals]------------------------------------!
  integer :: v, e, s
!==============================================================================!

  !-----------------------------!
  !   Initialize all elements   !
  !-----------------------------!
  do e = 1, MAX_SURFACE_ELEMENTS
    Front % elem(e) % nne  = 0
    Front % elem(e) % nv   = 0
    Front % elem(e) % ns   = 0
    Front % elem(e) % v(:) = 0
    Front % elem(e) % e(:) = 0
    Front % elem(e) % s(:) = 0
    Front % elem(e) % cell = 0
    Front % elem(e) % face = 0
  end do
  Front % n_elems = 0

  !-----------------------------!
  !   Initialize all vertices   !
  !-----------------------------!
  do v = 1, MAX_SURFACE_VERTICES

    Front % vert(v) % nne = 0

    ! Set initial velocity to zero
    Front % vert(v) % u_n = 0.0
    Front % vert(v) % v_n = 0.0
    Front % vert(v) % w_n = 0.0

    ! Set initial coordinates to zero
    Front % vert(v) % x_n = 0.0
    Front % vert(v) % y_n = 0.0
    Front % vert(v) % z_n = 0.0

    Front % vert(v) % x_o = 0.0
    Front % vert(v) % y_o = 0.0
    Front % vert(v) % z_o = 0.0

    ! Set initial cell, node and boundary cell to zero
    Front % vert(v) % cell     = 0
    Front % vert(v) % node     = 0
    Front % vert(v) % bnd_cell = 0
    Front % vert(v) % bnd_face = 0

    ! Assume vertex is in the domain
    ! (A smarter way could be worked out, depending ...
    ! ... on the result of the call to Find_Nearest_Cell)
    Front % vert(v) % escaped   = .false.

    ! Is vertex in this processor?
    Front % vert(v) % proc = 0
    Front % vert(v) % buff = 0

  end do
  Front % n_verts = 0

  !--------------------------!
  !   Initialize all sides   !
  !--------------------------!
  do s = 1, MAX_SURFACE_ELEMENTS * 3
    Front % side(s) % a        = 0
    Front % side(s) % b        = 0
    Front % side(s) % c        = 0
    Front % side(s) % d        = 0
    Front % side(s) % ei       = 0
    Front % side(s) % ea       = 0
    Front % side(s) % eb       = 0
    Front % side(s) % length   = 0.0
    Front % side(s) % boundary = .false.
  end do
  Front % n_sides = 0

  end subroutine
