!==============================================================================!
  subroutine Adjust_P_Drops(Flow, Grid)
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
!     acc = (u_bulk_o - u_bulk) / dt            [m/s / s = m/s^2]              !
!                                                                              !
!   * In a periodic channel, force is related to pressure drop as:             !
!                                                                              !
!     p_drop * volume = force                   [kg/m^2/s^2 * m^3 = N]         !
!                                                                              !
!   * Or:                                                                      !
!                                                                              !
!     p_drop = force / volume =                 [N/m^3        = kg/s^2/m^2]    !
!            = mass / volume * acc              [kg/m^3*m/s^2 = kg/s^2/m^2]    !
!            = rho * (u_bulk_o - u_bulk) / dt   [kg/m^3*m/s/s = kg/s^2/m^2]    !
!                                                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow  !! parent class for the flow
  type(Grid_Type)           :: Grid  !! grid on which the flow is defined
!-----------------------------------[Locals]-----------------------------------!
  type(Bulk_Type), pointer :: bulk
  real                     :: rho
!==============================================================================!

  ! Take aliases
  bulk => Flow % bulk

  ! Work out the mean density of the fluid body
  rho = Flow % Volume_Average(Grid, Flow % density)

  Assert(1.0 / Flow % dt > TINY)

  ! Once density is computed, continue with necessary pressure drops

  if(abs(bulk % u_o) >= TINY ) then
    bulk % p_drop_x = rho * (bulk % u_o - bulk % u) / Flow % dt
  end if
  if(abs(bulk % v_o) >= TINY ) then
    bulk % p_drop_y = rho * (bulk % v_o - bulk % v) / Flow % dt
  end if
  if(abs(bulk % w_o) >= TINY ) then
    bulk % p_drop_x = rho * (bulk % w_o - bulk % w) / Flow % dt
  end if

  end subroutine
