!==============================================================================!
  module Bulk_Mod
!------------------------------------------------------------------------------!
!   Volume fluxes, bulk velocities and pressure drops (for each material)      !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !   Bulk type   !
  !---------------!
  type Bulk_Type

    ! Cross-sectional areas defined by the monitoring planes
    real :: area_x  ! [m^2]
    real :: area_y  ! [m^2]
    real :: area_z  ! [m^2]

    ! Inflow and outflow areas
    real :: area_in   ! [m^2]
    real :: area_out  ! [m^2]

    ! Total inflow and outflow volume flux
    real :: vol_in   ! [m^3/s]
    real :: vol_out  ! [m^3/s]
    real :: vol_src  ! [m^3/s]

    ! Volume flow rates in x-, y-, z- direction
    real :: flux_x  ! [m^3/s]
    real :: flux_y  ! [m^3/s]
    real :: flux_z  ! [m^3/s]

    ! Desired Volume flow rates in x-, y-, z- direction (set in control file)
    real :: flux_x_o  ! [m^3/s]
    real :: flux_y_o  ! [m^3/s]
    real :: flux_z_o  ! [m^3/s]

    ! Pressure drop to achieve desired volume flux
    real :: p_drop_x  ! [N/m^3] = [kg/m^2/s^2]
    real :: p_drop_y  ! [N/m^3] = [kg/m^2/s^2]
    real :: p_drop_z  ! [N/m^3] = [kg/m^2/s^2]

    ! Bulk velocities
    real :: u  ! [m/s]
    real :: v  ! [m/s]
    real :: w  ! [m/s]

    ! Monitoring plane coordinates used to compute bulk fluxes
    real :: xp
    real :: yp
    real :: zp

  end type

  contains

#   include "Bulk_Mod/Adjust_P_Drops.f90"
#   include "Bulk_Mod/Monitoring_Planes_Areas.f90"
#   include "Bulk_Mod/Print_Areas.f90"

  end module
