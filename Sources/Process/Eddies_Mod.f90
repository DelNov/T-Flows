!==============================================================================!
  module Eddies_Mod
!------------------------------------------------------------------------------!
!>  The Eddies_Mod module is designed to manage the creation and manipulation
!>  of synthetic eddies, which are used in scale-resolving simulations like
!>  LES, DES and hybrid RANS-LES models. These simulations require a detailed
!>  representation of turbulent flow structures, and synthetic eddies provide
!>  a method to introduce these structures at the simulation's inlet boundary.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Math_Mod
  use Field_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !--------------------------------!
  !   A few important parameters   !
  !--------------------------------!
  integer, parameter :: MAX_TURBULENT_PLANES = 64

  !---------------!
  !   Eddy type   !
  !---------------!
  !> This type definition is used to describe a single eddy's characteristics.
  type Eddy_Type
    real :: x, y, z  !! eddy's position
    real :: u, v, w  !! eddy's velocity
    real :: length   !! eddy's length
    real :: radius   !! eddy's radius
    real :: sense    !! eddy's sense of rotation
  end type

  !-----------------!
  !   Eddies type   !
  !-----------------!
  !> This type is a collection (array) of individual eddies (Eddy_Type).
  !> It is associated with a specific grid (pnt_grid) and flow field
  !> (pnt_flow) and is defined for a particular boundary condition (bc_name).
  type Eddies_Type

    type(Grid_Type),  pointer :: pnt_grid  !! grid for which it is defined
    type(Field_Type), pointer :: pnt_flow  !! flow for which it is defined
    character(SL)             :: bc_name   !! boundary condition name

    integer                      :: n_eddies  !! number of eddies
    type(Eddy_Type), allocatable :: eddy(:)   !! storage for eddies

    real :: max_radius  !! maximum radius for an eddy
    real :: intensity   !! intensity of the eddies
    real :: x_plane
    real :: y_plane
    real :: z_plane

    ! Data specific to boundary condition plane eddies are defined
    integer :: n_bnd_cells      !! number of boundary cells in this processor
    integer :: n_bnd_cells_glo  !! number of boundary cells over all processors

    ! Coordinates, wall distance and velocity components
    ! of boundary cells over all processors
    ! (coordinates are also potential eddy centers)
    real, allocatable :: bnd_xc(:)
    real, allocatable :: bnd_yc(:)
    real, allocatable :: bnd_zc(:)
    real, allocatable :: bnd_wd(:)
    real, allocatable :: bnd_u (:)
    real, allocatable :: bnd_v (:)
    real, allocatable :: bnd_w (:)

    contains
      procedure          :: Advance_Eddies
      procedure          :: Create_Eddies
      procedure, private :: Gather_Bnd_Cells
      procedure, private :: Place_Eddy
      procedure          :: Superimpose_Eddies
  end type

  !--------------------------------------!
  !   Synthtetic turbulence plane type   !
  !--------------------------------------!
  !> This type represents a turbulent plane that can consist of multiple
  !> eddies. It includes the number of planes (n_planes) and an array of
  !> Eddies_Type (one for each turbulent plane).
  type Turb_Plane_Type
    integer           :: n_planes
    type(Eddies_Type) :: Plane(MAX_TURBULENT_PLANES)
  end type

  contains

#   include "Eddies_Mod/Advance.f90"
#   include "Eddies_Mod/Create.f90"
#   include "Eddies_Mod/Gather_Bnd_Cells.f90"
#   include "Eddies_Mod/Place_Eddy.f90"
#   include "Eddies_Mod/Superimpose.f90"

  end module
