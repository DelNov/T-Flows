!==============================================================================!
  module Line_Mod
!------------------------------------------------------------------------------!
!>  This module is responsible for representing the concept of a discrete line
!>  (1D grid) within the context of mesh generation in the sub-program Generate.
!>  Together with Point_Mod (Point_Type) and Block_Mod (Block_Type) it helps to
!>  define a computational domain in Generate's Domain_Mod.
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !   Line type   !
  !---------------!
  !> Encapsulates the properties and characteristics
  !> of a discrete line (1D grid) within the grid.
  type Line_Type

    integer           :: points(2)   !! strat and end points defining the line
    real, allocatable :: x(:)        !! x coordinates of points along the line
    real, allocatable :: y(:)        !! y coordinates of points along the line
    real, allocatable :: z(:)        !! z coordinates of points along the line
    real              :: weight      !! line weights for node clustering
    integer           :: resolution  !! line resolution (number of points)

  end type

  end module
