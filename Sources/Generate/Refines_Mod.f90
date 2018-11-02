!==============================================================================!
  module Refines_Mod
!------------------------------------------------------------------------------!
!   Type for refining a grid.                                                  !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Shapes of refinement
  integer, parameter :: ELIPSOID  = 3
  integer, parameter :: RECTANGLE = 4
  integer, parameter :: PLANE     = 5

  !------------------!
  !   Refines_Type   !
  !------------------!
  type Refines_Type

    integer              :: n_levels        ! number of refinement levels
    integer, allocatable :: cell_level(:)   ! cells' refinement level
    logical, allocatable :: cell_marked(:)  ! true if cell markered

    integer, allocatable :: n_regions(:)    ! number of refin. regions
    real,    allocatable :: region(:,:,:)   ! levels, regions

  end type

  !---------------------------!
  !   Member-like functions   !
  !---------------------------!
  contains

  include 'Refines_Mod/Allocate_Cells.f90'
  include 'Refines_Mod/Allocate_Levels.f90'

  end module
