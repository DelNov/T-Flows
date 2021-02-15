!==============================================================================!
  subroutine Front_Mod_Initialize(front)
!------------------------------------------------------------------------------!
!   Surface genesis                                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Front_Type), target :: front
!----------------------------------[Locals]------------------------------------!
  integer :: v, e, s
!==============================================================================!

  !-----------------------------!
  !   Initialize all elements   !
  !-----------------------------!
  do e = 1, MAX_SURFACE_ELEMENTS
    front % elem(e) % nne  = 0
    front % elem(e) % nv   = 0
    front % elem(e) % ns   = 0
    front % elem(e) % v(:) = 0
    front % elem(e) % e(:) = 0
    front % elem(e) % s(:) = 0
  end do
  front % n_elems = 0

  !-----------------------------!
  !   Initialize all vertices   !
  !-----------------------------!
  do v = 1, MAX_SURFACE_VERTICES

    front % vert(v) % nne = 0

    ! Set initial velocity to zero
    front % vert(v) % u_n = 0.0
    front % vert(v) % v_n = 0.0
    front % vert(v) % w_n = 0.0

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
