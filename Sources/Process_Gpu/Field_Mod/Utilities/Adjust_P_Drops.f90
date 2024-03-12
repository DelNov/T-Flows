!==============================================================================!
  subroutine Adjust_P_Drops(Flow)
!------------------------------------------------------------------------------!
!>  Adjusts pressure drops in the Bulk_Mod module to maintain constant volume
!>  flow rate in the computational domain. This subroutine is based on the
!>  Newton's first law and fluid dynamic principles, ensuring that the
!>  pressure drops are recalculated in response to changes in flow rates over
!>  time.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Recalculate the pressure drop to keep the constant volume flux           !
!                                                                              !
!   * First Newtons law:                                                       !
!                                                                              !
!     force = mass * acc                        [kg * m/s^2 = N]               !
!                                                                              !
!   * Necessary acceleration (to achieve the desired volume flux) is then:     !
!                                                                              !
!     acc = (v_flux_o - v_flux) / dt / area     [m^3/s / s / m^2 = m/s^2]      !
!                                                                              !
!   * In a periodic channel, force is related to pressure drop as:             !
!                                                                              !
!     p_drop * volume = force                   [kg/m^2/s^2 * m^3 = N]         !
!                                                                              !
!   * Or:                                                                      !
!                                                                              !
!     p_drop = force / volume =                                                !
!            = mass / volume * acc                                             !
!            = rho * (v_flux_o - v_flux) / dt / area                           !
!                                                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow  !! parent class
!-----------------------------------[Locals]-----------------------------------!
  type(Bulk_Type), pointer :: bulk
  type(Grid_Type), pointer :: Grid
  real                     :: rho
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  bulk => Flow % bulk

  ! Work out the mean density of the fluid body
  rho = Flow % density

  ! Once density is computed, continue with necessary pressure drops
  if( abs(bulk % flux_x_o) >= TINY ) then
    bulk % p_drop_x = rho * (bulk % flux_x_o - bulk % flux_x)  &
                    / (Flow % dt * bulk % area_x + TINY)
  end if
  if( abs(bulk % flux_y_o) >= TINY ) then
    bulk % p_drop_y = rho * (bulk % flux_y_o - bulk % flux_y)  &
                    / (Flow % dt * bulk % area_y + TINY)
  end if
  if( abs(bulk % flux_z_o) >= TINY ) then
    bulk % p_drop_z = rho * (bulk % flux_z_o - bulk % flux_z)  &
                    / (Flow % dt * bulk % area_z + TINY)
  end if

  end subroutine
