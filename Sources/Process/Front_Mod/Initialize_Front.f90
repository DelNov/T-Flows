!==============================================================================!
  subroutine Initialize_Front(Front)
!------------------------------------------------------------------------------!
!>  This subroutine is designed for initializing the front structure in a
!>  computational simulation. It focuses on setting up all the elements,
!>  vertices, and sides of the front. It iterates through each element, vertex,
!>  and side, setting their attributes to the default values.  The rimary role
!>  of this subroutine is to ensure that the front object starts from a clean,
!>  default state, ready for further processing and manipulation within the
!>  simulation environment.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type), target :: Front  !! parent class
!----------------------------------[Locals]------------------------------------!
  integer :: v, e, s
!==============================================================================!

  !-----------------------------!
  !   Initialize all elements   !
  !-----------------------------!
  do e = 1, size(Front % Elem, 1)
    call Front % Elem(e) % Initialize_Elem()
  end do
  Front % n_elems = 0

  !-----------------------------!
  !   Initialize all vertices   !
  !-----------------------------!
  do v = 1, size(Front % Vert, 1)
    call Front % Vert(v) % Initialize_Vert(Front % pnt_grid)
    Front % Vert(v) % trapped = .true.
  end do
  Front % n_verts = 0

  !--------------------------!
  !   Initialize all sides   !
  !--------------------------!
  do s = 1, size(Front % side, 1)
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
