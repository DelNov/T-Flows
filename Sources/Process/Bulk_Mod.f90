!==============================================================================!
  module Bulk_Mod
!------------------------------------------------------------------------------!
!   Mass fluxes, bulk velocities and pressure drops (for each material)        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Grid_Mod
  use Comm_Mod
  use Bnd_Cond_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !   Bulk type   !
  !---------------!
  type Bulk_Type

    real :: area_x
    real :: area_y
    real :: area_z

    real :: mass_in
    real :: mass_out

    ! Fluxes in x, y and z direction
    real :: flux_x
    real :: flux_y
    real :: flux_z

    ! Old fluxes in x, y and z direction
    real :: flux_x_o
    real :: flux_y_o
    real :: flux_z_o

    real :: p_drop_x
    real :: p_drop_y
    real :: p_drop_z

    real :: u
    real :: v
    real :: w

    ! Monitoring plane coordinates
    real :: xp
    real :: yp
    real :: zp

  end type

  contains

  include 'Bulk_Mod/Adjust_P_Drops.f90'
  include 'Bulk_Mod/Calculate_Fluxes.f90'
  include 'Bulk_Mod/Monitoring_Planes_Areas.f90'
  include 'Bulk_Mod/Print_Areas.f90'

  end module
