!==============================================================================!
  subroutine Initialize_Surf(Surf)
!------------------------------------------------------------------------------!
!>  This subroutine is designed for initializing the surface object in a
!>  computational simulation. It focuses on setting up all the elements,
!>  vertices, and sides of the surface. It iterates through each element,
!>  vertex, and side, setting their attributes to the default values.  The
!>  primary role of this subroutine is to ensure that the front object starts
!>  from a clean, default state, ready for further processing and manipulation
!>  within the simulation environment.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf  !! parent class
!----------------------------------[Locals]------------------------------------!
  integer :: v, e, s
!==============================================================================!

  !-----------------------------!
  !   Initialize all elements   !
  !-----------------------------!
  do e = 1, size(Surf % Elem, 1)
    call Surf % Elem(e) % Initialize_Elem()
  end do
  Surf % n_elems = 0

  !-----------------------------!
  !   Initialize all vertices   !
  !-----------------------------!
  do v = 1, size(Surf % Vert, 1)
    call Surf % Vert(v) % Initialize_Vert(Surf % pnt_grid)
    Surf % Vert(v) % trapped = .true.
  end do
  Surf % n_verts = 0

  !--------------------------!
  !   Initialize all sides   !
  !--------------------------!
  do s = 1, size(Surf % side, 1)
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
