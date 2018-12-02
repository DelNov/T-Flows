!==============================================================================!
  subroutine Read_Physical(grid, restar)
!------------------------------------------------------------------------------!
!   Reads details about physial models.                                        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod
  use Rans_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  logical         :: restar
!==============================================================================!

  !-------------------------!
  !   Related to bouyancy   !
  !-------------------------!
  call Control_Mod_Gravitational_Vector(.true.)
  call Control_Mod_Buoyancy(.true.)
  call Control_Mod_Reference_Temperature(.true.)

  !---------------------------!
  !   Related to turbulence   !
  !---------------------------!
  call Control_Mod_Turbulence_Model(.true.)
  call Control_Mod_Turbulence_Model_Variant(.true.)
  call Control_Mod_Rough_Walls(.true.)
  call Control_Mod_Turbulent_Heat_Flux_Model(.true.)

  !-------------------------------------------------------------------------!
  !   Initialization of model constants depending on the turbulence model   !
  !-------------------------------------------------------------------------!
  if(turbulence_model .eq. K_EPS) then
    call Constants_K_Eps()
  end if

  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
    call Constants_Reynolds_Stress()
  end if

  if(turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
    call Constants_Hanjalic_Jakirlic()
  end if

  if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
     turbulence_model .eq. HYBRID_LES_RANS) then
    call Constants_K_Eps_Zeta_F()
  end if

  if(turbulence_model .eq. SPALART_ALLMARAS .or.  &
     turbulence_model .eq. DES_SPALART) then
    call Constants_Spalart_Allmaras()
  end if

  !------------------------------------!
  !   Pressure drops and mass fluxes   !
  !------------------------------------!
  if(.not. restar) then
    call Control_Mod_Pressure_Drops()
  end if
  if(.not. restar) then
    call Control_Mod_Mass_Flow_Rates()
  end if

  end subroutine
