!==============================================================================!
  subroutine Initialize_Vert(Vert, Grid)
!------------------------------------------------------------------------------!
!   Initializes point, as the name clearly implies                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vert_Type)        :: Vert
  type(Grid_Type), target :: Grid
!==============================================================================!

  ! Call parent's constructor
  call Vert % Initialize_Point(Grid)

  ! Old vertex's coordinates; new is in the parent
  Vert % x_o = 0.0
  Vert % y_o = 0.0
  Vert % z_o = 0.0

  ! Vertex's new velocity (if ever needed)
  Vert % u_n = 0.0
  Vert % v_n = 0.0
  Vert % w_n = 0.0

  ! I forgot what these were, something for smoothing, I think
  Vert % sumx = 0.0
  Vert % sumy = 0.0
  Vert % sumz = 0.0

  Vert % nne      = 0        ! number of neighbouring elements
  Vert % nnv      = 0        ! number of neighbouring vertices
  Vert % boundary = .false.  ! is vertex on a boundary
  Vert % e(1:24)  = 0        ! list of elements around the vertex

  ! The closest cell, node, boundary cell and face
  Vert % bnd_face = 0

  ! Vertex departure from domain
  Vert % escaped = .false.

  ! Curvature at the vertex
  Vert % curv = 0.0

  end subroutine
