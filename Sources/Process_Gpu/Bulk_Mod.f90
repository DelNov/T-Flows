#include "../Shared/Browse.h90"

!==============================================================================!
  module Bulk_Mod
!------------------------------------------------------------------------------!
!>  The Bulk_Mod module in T-Flows' Process originated as a utility to maintain
!>  consistent volume flow rates and, associated with that, bulk velocities,
!>  and Reynolds numbers for LES studies of flows in channels and cube matrices.
!>  The module is not limited to these flows, and is effective for any flow
!>  with periodicity in the streamwise direction, providing a systematic
!>  approach to keep the desire volume flow rates or monitor flow rates as a
!>  function of pressure drops.  During a CFD simulation, Process even prints
!>  the current values of volume flow rates and pressure drops on the terminal.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !   Bulk type   !
  !---------------!
  !> Bulk_Type encapsulates data which facilitates computation of volume
  !> fluxes through a computational domain and adjust pressure drops if
  !> instructed to do so.
  type Bulk_Type

    ! Inflow and outflow areas
    real :: area_in   !! inflow area [m^2]
    real :: area_out  !! outflow area [m^2]

    ! Total inflow and outflow volume flow rate
    real :: vol_in    !! volume flow rate coming from inflows [m^3/s]
    real :: vol_out   !! volume flow rate exiting through outflows [m^3/s]
    real :: vol_src   !! internal volumetric source due to phase change [m^3/s]

    ! Pressure drop to achieve the target volume flow rates
    real :: p_drop_x  !! pressure drop in x direction [N/m^3] = [kg/m^2/s^2]
    real :: p_drop_y  !! pressure drop in y direction [N/m^3] = [kg/m^2/s^2]
    real :: p_drop_z  !! pressure drop in z direction [N/m^3] = [kg/m^2/s^2]

    ! Bulk velocities
    real :: u         !! x component of bulk velocity [m/s]
    real :: v         !! y component of bulk velocity [m/s]
    real :: w         !! z component of bulk velocity [m/s]

    ! Target bulk velocities
    real :: u_o       !! x component of target bulk velocity [m/s]
    real :: v_o       !! y component of target bulk velocity [m/s]
    real :: w_o       !! z component of target bulk velocity [m/s]

  end type

  end module
