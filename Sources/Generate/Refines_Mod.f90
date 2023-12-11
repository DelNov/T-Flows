#include "../Shared/Assert.h90"

!==============================================================================!
  module Refines_Mod
!------------------------------------------------------------------------------!
!>  The Refines_Mod module in the Generate program defines a structure and
!>  functionality for refining a computational grid.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Point_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Shapes of refinement
  integer, parameter :: ELIPSOID  = 3  !! elipsoidal refinement region
  integer, parameter :: RECTANGLE = 4  !! rectangular refinement region
  integer, parameter :: PLANE     = 5  !! refinement region defined by a plane

  ! Maximum number of shapes
  integer, parameter :: MAX_SHAPES = 1024  !! max number of refinement regions

  !----------------!
  !   Shape_Type   !
  !----------------!
  !> A type that encapsulates the shape of a refinement region.
  type Shape_Type
    integer          :: shape   !! shape tag (ELLIPSOID, RECTANGLE, PLANE)
    type(Point_Type) :: pnt(2)  !! first and second point defining the position
  end type

  !------------------!
  !   Refines_Type   !
  !------------------!
  !> Encapsulates data pertinent to cell refinement
  type Refines_Type

    integer              :: n_levels        !! number of refinement levels
    integer, allocatable :: cell_level(:)   !! cells' refinement level
    logical, allocatable :: cell_marked(:)  !! true if cell markered

    integer,          allocatable :: n_ranges(:)  !! number refinement regions
    type(Shape_Type), allocatable :: range(:,:)   !! over levels and regions

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
