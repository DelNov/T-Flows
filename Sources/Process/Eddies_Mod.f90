!==============================================================================!
  module Eddies_Mod
!------------------------------------------------------------------------------!
!   Module to define eddy type used in synthetic eddy generation               !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Math_Mod
  use Grid_Mod
  use Comm_Mod
  use Field_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

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
    character(len=80)         :: bc_name   ! boundary condition name

    integer                      :: n_eddies  ! number of eddies
    type(Eddy_Type), allocatable :: eddy(:)   ! storage for eddies

    real :: max_radius  ! maximum radius for an eddy
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

  contains

  include 'Eddies_Mod/Advance.f90'
  include 'Eddies_Mod/Allocate.f90'
  include 'Eddies_Mod/Gather_Bnd_Cells.f90'
  include 'Eddies_Mod/Place_Eddy.f90'
  include 'Eddies_Mod/Superimpose.f90'

  end module
