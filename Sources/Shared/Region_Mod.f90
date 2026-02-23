!==============================================================================!
  module Region_Mod
!------------------------------------------------------------------------------!
!>  The Region_Mod module in T-Flows is designed to manage boundary condition
!>  regions within the Grid_Type. It categorizes various types of boundary
!>  conditions, providing a structured approach to handle different boundary
!>  interactions in computational domains. The module is also used to store
!>  threads for OMP runs, although OMP is in an early stage of development.
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
  integer, parameter :: AMBIENT    = 10079
  integer, parameter :: INSIDE     = 10091
  integer, parameter :: BUFFER     = 10093
  integer, parameter :: PERIODIC_X = 10099
  integer, parameter :: PERIODIC_Y = 10103
  integer, parameter :: PERIODIC_Z = 10111
  integer, parameter :: UNDEFINED  = 10133

  !-----------------!
  !   Region type   !
  !-----------------!
  !> Region_Type hold ariables for storing boundary condition information
  !> such as names and types of regions, regions at cells and faces and
  !> first and last cell or face for each region for more efficient browsing.
  type Region_Type

    ! Name of the boundary conditions specified in grid generation
    ! It ranges through number of boundary conditions (aka regions)
    character(SL), allocatable :: name(:)  !! name of the boundary condition

    ! Boundary types, ranging through all regions
    integer, allocatable :: type(:)  !! type of the boundary condition

    ! Boundary condition ranging through boundary cells.
    ! Values start from one, zero is internal cell
    ! (Follows nomenclature from "../Shared/Comm_Mod_Par.f90")
    integer, allocatable :: at_cell(:)  !! region at cell
    integer, allocatable :: at_face(:)  !! region at face (for periodicity)
    integer, allocatable :: f_cell(:)   !! first (bnd) cell for region
    integer, allocatable :: l_cell(:)   !! last (bnd) cell for region
    integer, allocatable :: f_face(:)   !! first (bnd) cell for region
    integer, allocatable :: l_face(:)   !! last (bnd) cell for region

  end type

  end module
