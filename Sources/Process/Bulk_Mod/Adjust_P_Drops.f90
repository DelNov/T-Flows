!==============================================================================!
  subroutine Bulk_Mod_Adjust_P_Drops(bulk, dt)
!------------------------------------------------------------------------------!
!   Recalculate the pressure drop to keep the constant mass flux               !
!                                                                              !
!   First Newtons law:                                                         !
!                                                                              !
!   F = m * a                                                                  !
!                                                                              !
!   where:                                                                     !
!                                                                              !
!   a = dv / dt = dFlux / dt * 1 / (A * rho)                                   !
!   m = rho * v                                                                !
!   F = Pdrop * l * A = Pdrop * v                                              !
!                                                                              !
!   finally:                                                                   !
!                                                                              !
!   Pdrop * v = rho * v * dFlux / dt * 1 / (A * rho)                           !
!                                                                              !
!   after cancelling: v and rho, it yields:                                    !
!                                                                              !
!   Pdrop = dFlux/dt/A                                                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Bulk_Type) :: bulk
  real            :: dt
!==============================================================================!

  if( abs(bulk % flux_x_o) >= TINY ) then
    bulk % p_drop_x = (bulk % flux_x_o - bulk % flux_x) &
                       / (dt * bulk % area_x + TINY)
  end if
  if( abs(bulk % flux_y_o) >= TINY ) then
    bulk % p_drop_y = (bulk % flux_y_o - bulk % flux_y) &
                       / (dt * bulk % area_y + TINY)
  end if
  if( abs(bulk % flux_z_o) >= TINY ) then
    bulk % p_drop_z = (bulk % flux_z_o - bulk % flux_z) &
                       / (dt * bulk % area_z + TINY)
  end if

  end subroutine
