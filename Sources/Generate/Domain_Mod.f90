!==============================================================================!
  module Domain_Mod
!------------------------------------------------------------------------------!
!   Domain as the one used in "Generator"                                      !
!------------------------------------------------------------------------------!
  use Gen_Mod                      ! a relict from the past
  use Point_Mod, only: Point_Type
  use Line_Mod,  only: Line_Type
  use Block_Mod, only: Block_Type
  use Range_Mod, only: Range_Type  ! to store boundary conditions
  use Grid_Mod
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
    integer :: n_ranges
    integer :: n_smooths

    type(Point_Type), allocatable :: points(:)
    type(Block_Type), allocatable :: blocks(:)
    type(Line_Type),  allocatable :: lines(:)
    type(Range_Type), allocatable :: ranges(:)

  end type

  !---------------------------!
  !   Member-like functions   !
  !---------------------------!
  contains

#   include "Domain_Mod/Allocate_Points.f90"
#   include "Domain_Mod/Allocate_Blocks.f90"
#   include "Domain_Mod/Allocate_Lines.f90"
#   include "Domain_Mod/Allocate_Ranges.f90"
#   include "Domain_Mod/Calculate_Node_Coordinates.f90"
#   include "Domain_Mod/Connect_Blocks.f90"
#   include "Domain_Mod/Connect_Periodicity.f90"
#   include "Domain_Mod/Distribute_Nodes.f90"
#   include "Domain_Mod/Distribute_Ranges.f90"
#   include "Domain_Mod/Find_Line.f90"
#   include "Domain_Mod/Find_Surface.f90"
#   include "Domain_Mod/Is_Line_In_Block.f90"
#   include "Domain_Mod/Laplace.f90"

  end module
