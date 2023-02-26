!==============================================================================!
  module Region_Mod
!------------------------------------------------------------------------------!
!   This is used to store regions (boundary conditions) within a Grid_Type     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------------------------------------------------!
  !   Constants for identification of boundary conditions   !
  !---------------------------------------------------------!
  integer, parameter :: INFLOW     = 10007
  integer, parameter :: WALL       = 10009
  integer, parameter :: OUTFLOW    = 10037
  integer, parameter :: SYMMETRY   = 10039
  integer, parameter :: CONVECT    = 10061
  integer, parameter :: WALLFL     = 10067
  integer, parameter :: PRESSURE   = 10069
  integer, parameter :: INSIDE     = 10079
  integer, parameter :: BUFFER     = 10091
  integer, parameter :: PERIODIC_X = 10093
  integer, parameter :: PERIODIC_Y = 10099
  integer, parameter :: PERIODIC_Z = 10103
  integer, parameter :: UNDEFINED  = 10111

  !-----------------!
  !   Region type   !
  !-----------------!
  type Region_Type

    ! Name of the boundary conditions specified in grid generation
    ! It ranges through number of boundary conditions (aka regions)
    character(SL), allocatable :: name(:)

    ! Boundary types, ranging through all regions
    integer, allocatable :: type(:)

    ! Boundary condition ranging through boundary cells.
    ! Values start from one, zero is internal cell
    ! (Follows nomenclature from "../Shared/Comm_Mod_Par.f90")
    integer, allocatable :: at_cell(:)  ! region at cell
    integer, allocatable :: f_cell(:)   ! first (bnd) cell for region
    integer, allocatable :: l_cell(:)   ! last (bnd) cell for region
    integer, allocatable :: f_face(:)   ! first (bnd) cell for region
    integer, allocatable :: l_face(:)   ! last (bnd) cell for region

  end type

  end module
