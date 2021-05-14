!==============================================================================!
  subroutine Read_Control_Physical_Properties(flow, mult, swarm)
!------------------------------------------------------------------------------!
!   Reads physical properties from control file.                               !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
  use Field_Mod
  use Swarm_Mod
  use Control_Mod
  use Multiphase_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type)      :: flow
  type(Multiphase_Type) :: mult
  type(Swarm_Type)      :: swarm
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  real                     :: dens_const, visc_const
  real                     :: capa_const, cond_const
!==============================================================================!

  ! Take alias
  grid => flow % pnt_grid

  ! Read constant (defualt) values
  call Control_Mod_Dynamic_Viscosity   (visc_const)
  call Control_Mod_Mass_Density        (dens_const)
  call Control_Mod_Heat_Capacity       (capa_const)
  call Control_Mod_Thermal_Conductivity(cond_const)
  call Control_Mod_Scalars_Diffusivity (flow % diffusivity)

  if(mult % model .eq. VOLUME_OF_FLUID) then
    call Control_Mod_Phase_Densities       (mult % phase_dens)
    call Control_Mod_Phase_Viscosities     (mult % phase_visc)
    call Control_Mod_Phase_Capacities      (mult % phase_capa)
    call Control_Mod_Phase_Conductivities  (mult % phase_cond)
    call Control_Mod_Surface_Tension       (mult % surface_tension)
    call Control_Mod_Latent_Heat           (mult % latent_heat)
    call Control_Mod_Saturation_Temperature(mult % t_sat)
  else
    flow % density     (:) = dens_const
    flow % viscosity   (:) = visc_const
    flow % capacity    (:) = capa_const
    flow % conductivity(:) = cond_const
  end if

  end subroutine
