!==============================================================================!
  module Refines_Mod
!------------------------------------------------------------------------------!
!   Type for refining a grid.                                                  !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Point_Mod
  use Gen_Mod    ! artifact of the past
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Shapes of refinement
  integer, parameter :: ELIPSOID  = 3
  integer, parameter :: RECTANGLE = 4
  integer, parameter :: PLANE     = 5

  ! Maximum number of shapes
  integer, parameter :: MAX_SHAPES = 1024

  !----------------!
  !   Shape_Type   !
  !----------------!
  type Shape_Type
    integer          :: shape   ! shape tag (ELLIPSOID, RECTANGLE, PLANE)
    type(Point_Type) :: pnt(2)  ! first and second point defining the position
  end type

  !------------------!
  !   Refines_Type   !
  !------------------!
  type Refines_Type

    integer              :: n_levels        ! number of refinement levels
    integer, allocatable :: cell_level(:)   ! cells' refinement level
    logical, allocatable :: cell_marked(:)  ! true if cell markered

    integer,          allocatable :: n_regions(:)  ! number refin. regions
    type(Shape_Type), allocatable :: region(:,:)   ! over levels & regions

  end type

  !---------------------------!
  !   Member-like functions   !
  !---------------------------!
  contains

#   include "Refines_Mod/Allocate_Cells.f90"
#   include "Refines_Mod/Allocate_Levels.f90"
#   include "Refines_Mod/Connectivity.f90"
#   include "Refines_Mod/Mark_Cells.f90"
#   include "Refines_Mod/Refine_Marked_Cells.f90"
#   include "Refines_Mod/Which_Node.f90"

  end module
