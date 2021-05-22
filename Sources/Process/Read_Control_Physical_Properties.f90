!==============================================================================!
  subroutine Read_Control_Physical_Properties(Flow, Vof, swarm)
!------------------------------------------------------------------------------!
!   Reads physical properties from control file.                               !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
  use Field_Mod
  use Swarm_Mod
  use Control_Mod
  use Vof_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type) :: Flow
  type(Vof_Type)   :: Vof
  type(Swarm_Type) :: swarm
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  real                     :: dens_const, visc_const
  real                     :: capa_const, cond_const
!==============================================================================!

  ! Take alias
  grid => Flow % pnt_grid

  ! Read constant (defualt) values
  call Control_Mod_Dynamic_Viscosity   (visc_const)
  call Control_Mod_Mass_Density        (dens_const)
  call Control_Mod_Heat_Capacity       (capa_const)
  call Control_Mod_Thermal_Conductivity(cond_const)
  call Control_Mod_Scalars_Diffusivity (Flow % diffusivity)

  if(Flow % with_interface) then
    call Control_Mod_Phase_Densities       (Vof % phase_dens)
    call Control_Mod_Phase_Viscosities     (Vof % phase_visc)
    call Control_Mod_Phase_Capacities      (Vof % phase_capa)
    call Control_Mod_Phase_Conductivities  (Vof % phase_cond)
    call Control_Mod_Surface_Tension       (Vof % surface_tension)
    call Control_Mod_Latent_Heat           (Vof % latent_heat)
    call Control_Mod_Saturation_Temperature(Vof % t_sat)
  else
    Flow % density     (:) = dens_const
    Flow % viscosity   (:) = visc_const
    Flow % capacity    (:) = capa_const
    Flow % conductivity(:) = cond_const
  end if

  end subroutine
