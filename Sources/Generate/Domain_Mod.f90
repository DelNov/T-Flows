!==============================================================================!
  module Domain_Mod
!------------------------------------------------------------------------------!
!   Domain as the one used in "Generator"                                      !
!------------------------------------------------------------------------------!
  use Point_Mod,  only: Point_Type
  use Line_Mod,   only: Line_Type
  use Block_Mod,  only: Block_Type
  use Region_Mod, only: Region_Type  ! to store boundary conditions
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-----------------!
  !   Domain type   !
  !-----------------!
  type Domain_Type

    integer :: n_points
    integer :: n_blocks
    integer :: n_lines
    integer :: n_regions
    integer :: n_smooths

    type(Point_Type),  allocatable :: points(:)
    type(Block_Type),  allocatable :: blocks(:)
    type(Line_Type),   allocatable :: lines(:)
    type(Region_Type), allocatable :: regions(:)

  end type

  !---------------------------!
  !   Member-like functions   !
  !---------------------------!
  contains

  include 'Domain_Mod/Allocate_Points.f90'
  include 'Domain_Mod/Allocate_Blocks.f90'
  include 'Domain_Mod/Allocate_Lines.f90'
  include 'Domain_Mod/Allocate_Regions.f90'

  end module
