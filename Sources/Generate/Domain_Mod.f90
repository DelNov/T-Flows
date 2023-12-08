!==============================================================================!
  module Domain_Mod
!------------------------------------------------------------------------------!
!>  The Domain_Mod module, as part of the Generate program is an essential part
!>  of the code. It integrates several specialized modules (Point_Mod, Line_Mod,
!>  Block_Mod and Range_Mod) to define a computational grid.
!------------------------------------------------------------------------------!
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
  !> Central to the Domain_Mod module, defining the computational domain.
  type Domain_Type

    integer :: n_points   !! number of points specified in the .dom file,
                          !! hence number of points defining the domain
    integer :: n_blocks   !! number of hexahedral blocks defining the domain
    integer :: n_lines    !! number of lines defined in the .dom file
    integer :: n_ranges   !! number of refinement regions
    integer :: n_smooths  !! number of smoothing regions

    integer              :: n_periodic_cond       !! number of periodic bounds.
    integer, allocatable ::   periodic_cond(:,:)  !! store periodic conditions

    type(Point_Type), allocatable :: points(:)  !! array to store points
    type(Block_Type), allocatable :: blocks(:)  !! array to store blocks
    type(Line_Type),  allocatable :: lines(:)   !! array to store lines
    type(Range_Type), allocatable :: ranges(:)  !! array to store ranges

    contains
      procedure :: Allocate_Points
      procedure :: Allocate_Blocks
      procedure :: Allocate_Lines
      procedure :: Allocate_Ranges
      procedure :: Calculate_Node_Coordinates
      procedure :: Connect_Blocks
      procedure :: Connect_Periodicity
      procedure :: Distribute_Nodes
      procedure :: Distribute_Ranges
      procedure :: Find_Line
      procedure :: Find_Surface
      procedure :: Is_Line_In_Block
      procedure :: Laplace

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
