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
    call Front % Elem(e) % Initialize_Elem()
  end do
  Front % n_elems = 0

  !-----------------------------!
  !   Initialize all vertices   !
  !-----------------------------!
  do v = 1, MAX_SURFACE_VERTICES
    call Front % Vert(v) % Initialize_Vert(Front % pnt_grid)
    Front % Vert(v) % trapped = .true.
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
