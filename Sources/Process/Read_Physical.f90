!==============================================================================!
  subroutine Read_Physical(grid, restar)
!------------------------------------------------------------------------------!
!   Reads details about physial models.                                        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Flow_Mod
  use Rans_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  logical         :: restar
!-----------------------------------[Locals]-----------------------------------!
  integer :: m
!==============================================================================!

  call Control_Mod_Gravitational_Vector(grav_x, grav_y, grav_z)

  call Control_Mod_Turbulence_Model(.true.)

  call Control_Mod_Turbulence_Model_Variant(.true.)

  call Control_Mod_Rough_Walls(.true.)

  if(turbulence_model .eq. K_EPS) then
    call Constants_K_Eps()
  end if

  if(turbulence_model .eq. REYNOLDS_STRESS) then
    call Constants_Reynolds_Stress()
  end if

  if(turbulence_model .eq. HANJALIC_JAKIRLIC) then
    call Constants_Hanjalic_Jakirlic()
  end if

  if(turbulence_model .eq. K_EPS_ZETA_F) then
    call Constants_K_Eps_Zeta_F()
  end if

  if(turbulence_model .eq. SPALART_ALLMARAS .or.  &
     turbulence_model .eq. DES_SPALART) then
    call Constants_Spalart_Allmaras()
  end if

  ! Pressure drops
  if(.not. restar) then
    call Control_Mod_Pressure_Drops(bulk % p_drop_x,  &
                                    bulk % p_drop_y,  &
                                    bulk % p_drop_z)
  end if

  ! Mass fluxes
  if(.not. restar) then
    call Control_Mod_Mass_Flow_Rates(bulk % flux_x_o,  &
                                     bulk % flux_y_o,  &
                                     bulk % flux_z_o)
  end if

  end subroutine
