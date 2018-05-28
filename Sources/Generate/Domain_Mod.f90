!==============================================================================!
  module Domain_Mod
!------------------------------------------------------------------------------!
!   Domain as the one used in "Generator"                                      !
!------------------------------------------------------------------------------!
  use Point_Mod
  use Line_Mod
  use Block_Mod
  use Region_Mod  ! to store boundary conditions or materials
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
