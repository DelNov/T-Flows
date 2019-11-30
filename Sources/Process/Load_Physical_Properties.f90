!==============================================================================!
  subroutine Load_Physical_Properties(flow)
!------------------------------------------------------------------------------!
!   Reads physical properties from control file.                               !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
  use Field_Mod
  use Control_Mod
  use Multiphase_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type) :: flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  real                     :: dens_const, visc_const
!==============================================================================!

  ! Take alias
  grid => flow % pnt_grid

  ! Allocate properties
  allocate(flow % density  (-grid % n_bnd_cells:grid % n_cells))
  allocate(flow % viscosity(-grid % n_bnd_cells:grid % n_cells))
  allocate(flow % density_f( grid % n_faces))


  call Control_Mod_Dynamic_Viscosity   (visc_const)
  call Control_Mod_Mass_Density        (dens_const)
  call Control_Mod_Heat_Capacity       (flow % capacity)
  call Control_Mod_Thermal_Conductivity(flow % conductivity)
  call Control_Mod_Species_Diffusivity (flow % diffusivity)

  if(multiphase_model .eq. VOLUME_OF_FLUID) then
    call Control_Mod_Phase_Densities  (phase_dens)
    call Control_Mod_Phase_Viscosities(phase_visc)
    call Control_Mod_Surface_Tension  (surface_tension)
  else
    flow % density(:)   = dens_const
    flow % viscosity(:) = visc_const
    flow % density_f(:) = dens_const
  end if

  end subroutine
