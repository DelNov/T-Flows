!==============================================================================!
  module Point_Mod
!------------------------------------------------------------------------------!
!   The parent of all point-based modules such as vertices and particles       !
!   Not to be confused with points in a numerical grid, those are called       !
!   nodes, and do not originate from this class.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !----------------!
  !   Point type   !
  !----------------!
  type Point_Type

    ! Grid in which the point is defined
    type(Grid_Type),  pointer :: pnt_grid

    ! Point's coordinates; new only, a point has no concept of history
    real :: x_n
    real :: y_n
    real :: z_n

    ! The closest cell, node, boundary cell (and face?)
    ! in the grid in which the point is defined
    integer :: cell
    integer :: node
    integer :: bnd_cell

    ! Point inside the subdomain
    integer :: proc
    integer :: buff

    contains
      procedure :: Find_Nearest_Cell
      procedure :: Find_Nearest_Node
      procedure :: Initialize_Point

  end type

  contains

#   include "Point_Mod/Find_Nearest_Cell.f90"
#   include "Point_Mod/Find_Nearest_Node.f90"
#   include "Point_Mod/Initialize_Point.f90"

  end module
