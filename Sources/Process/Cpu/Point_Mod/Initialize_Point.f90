!==============================================================================!
  subroutine Initialize_Point(Point, Grid)
!------------------------------------------------------------------------------!
!>  This subroutine is responsible for initializing a point within the context
!>  of a computational grid. It sets the point's basic properties to default
!>  values, ensuring a clean starting state for further operations.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Grid association: Establishes the relationship between the point and the !
!     grid within which it is defined.                                         !
!   * Coordinate initialization: Sets the initial coordinates of the point to  !
!     default values (typically zeroes).                                       !
!   * Proximity initialization: Initializes the closest cell, node, and        !
!     boundary cell indicators to default values.                              !
!   * Processor and buffer assignment: Initializes the processor and buffer    !
!     identifiers for the point to default values.                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Point_Type)       :: Point  !! point being initialized
  type(Grid_Type), target :: Grid   !! grid for which the point is defined
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
