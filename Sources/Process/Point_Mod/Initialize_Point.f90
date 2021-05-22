!==============================================================================!
  subroutine Initialize_Point(Point, Grid)
!------------------------------------------------------------------------------!
!   Initializes point, as the name clearly implies                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Point_Type)       :: Point
  type(Grid_Type), target :: Grid
!==============================================================================!

  ! Pointer to the grid for which the point is defined
  Point % pnt_grid => Grid

  ! Point's coordinates
  Point % x_n = 0.0
  Point % y_n = 0.0
  Point % z_n = 0.0

  ! Point's closes cell, node and boundary cell
  Point % cell     = 0
  Point % node     = 0
  Point % bnd_cell = 0

  ! The processor, and the buffer in which the point resides
  Point % proc = 0
  Point % buff = 0

  end subroutine
