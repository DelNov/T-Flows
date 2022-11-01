!==============================================================================!
  module Eddies_Mod
!------------------------------------------------------------------------------!
!   Module to define eddy type used in synthetic eddy generation               !
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
  !               !
  !   Eddy type   !
  !               !
  !---------------!
  type Eddy_Type
    real :: x, y, z  ! eddy's position
    real :: u, v, w  ! eddy's velocity
    real :: length   ! eddy's length
    real :: radius   ! eddy's radius
    real :: sense    ! eddy's sense of rotation
  end type

  !-----------------!
  !                 !
  !   Eddies type   !
  !                 !
  !-----------------!
  type Eddies_Type

    type(Grid_Type),  pointer :: pnt_grid  ! grid for which it is defined
    type(Field_Type), pointer :: pnt_flow  ! flow for which it is defined
    character(SL)             :: bc_name   ! boundary condition name

    integer                      :: n_eddies  ! number of eddies
    type(Eddy_Type), allocatable :: eddy(:)   ! storage for eddies

    real :: max_radius  ! maximum radius for an eddy
    real :: intensity   ! intensity of the eddies
    real :: x_plane
    real :: y_plane
    real :: z_plane

    ! Data specific to boundary condition plane eddies are defined
    integer :: n_bnd_cells      ! number of boundary cells in this processor
    integer :: n_bnd_cells_glo  ! number of boundary cells over all processors

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
  end type

  !--------------------------------------!
  !                                      !
  !   Synthtetic turbulence plane type   !
  !                                      !
  !--------------------------------------!
  type Turb_Plane_Type
    integer           :: n_planes
    type(Eddies_Type) :: plane(MAX_TURBULENT_PLANES)
  end type

  contains

#   include "Eddies_Mod/Advance.f90"
#   include "Eddies_Mod/Allocate.f90"
#   include "Eddies_Mod/Gather_Bnd_Cells.f90"
#   include "Eddies_Mod/Place_Eddy.f90"
#   include "Eddies_Mod/Superimpose.f90"

  end module
