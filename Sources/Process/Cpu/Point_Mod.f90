#include "../../Shared/Browse.h90"

!==============================================================================!
  module Point_Mod
!------------------------------------------------------------------------------!
!>  The Point_Mod module serves as the parent for all point-based modules in
!>  T-Flows, such as vertices and particles. It is distinct from numerical grid
!>  points (referred to as nodes) and is exclusively used within the Process
!>  sub-program. The module provides functionality to manage and interact with
!>  points in the computational domain, offering procedures to determine their
!>  nearest cells, nodes, and boundary cells.  A module with the same name
!>  exists in Generate, but that serves a different purpose from this one.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !----------------!
  !   Point type   !
  !----------------!
  !> Point_Type includes basic properties of points such as coordinates,
  !> the closest cell, node, boundary cell, and information about the
  !> point's position within a subdomain / buffer.
  type Point_Type

    ! Grid in which the point is defined
    type(Grid_Type),  pointer :: pnt_grid  !! pointer to computational grid

    ! Point's coordinates; new only, a point has no concept of history
    real :: x_n  !! point's x coordinate
    real :: y_n  !! point's y coordinate
    real :: z_n  !! point's z coordinate

    ! The closest cell, node, boundary cell (and face?)
    ! in the grid in which the point is defined
    integer :: cell      !! computational cell closest to the point
    integer :: node      !! computational node closest to the point
    integer :: bnd_cell  !! boundary cell closest to the point

    ! Point inside the subdomain
    integer :: proc  !! in which processor does the point reside
    integer :: buff  !! in which buffer cell does the point reside

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
