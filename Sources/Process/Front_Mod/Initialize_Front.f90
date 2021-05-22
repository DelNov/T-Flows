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
    Front % elem(e) % nx   = 0.0
    Front % elem(e) % ny   = 0.0
    Front % elem(e) % nz   = 0.0
    Front % elem(e) % xc   = 0.0
    Front % elem(e) % yc   = 0.0
    Front % elem(e) % zc   = 0.0
    Front % elem(e) % xe   = 0.0
    Front % elem(e) % ye   = 0.0
    Front % elem(e) % ze   = 0.0
    Front % elem(e) % sx   = 0.0
    Front % elem(e) % sy   = 0.0
    Front % elem(e) % sz   = 0.0
  end do
  Front % n_elems = 0

  !-----------------------------!
  !   Initialize all vertices   !
  !-----------------------------!
  do v = 1, MAX_SURFACE_VERTICES
    call Front % Vert(v) % Initialize_Vert(Front % pnt_grid)
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
