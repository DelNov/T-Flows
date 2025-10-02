!==============================================================================!
  subroutine Physical_Properties(Rc, Flow, Vof, Swarm)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to read the control file and set up the
!>  physical properties of the fluid for single and multi-phase flow
!>  simulations.  The subroutine focuses only on constant physical properties.
!>  Any variation in physical properties should be adressed in user functions.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Read_Controls_Type), intent(in) :: Rc     !! parent class
  type(Field_Type)                      :: Flow   !! flow object
  type(Vof_Type)                        :: Vof    !! VOF object
  type(Swarm_Type)                      :: Swarm  !! swarm object
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  real                     :: dens_const, visc_const
  real                     :: capa_const, cond_const
  real                     :: diff_const
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Rc)
  Unused(Swarm)
!==============================================================================!

  ! Give some sign
  O_Print '(a)', ' # Reading about physical properties'

  ! Take alias
  Grid => Flow % pnt_grid

  ! Read constant (defualt) values
  call Control % Dynamic_Viscosity(visc_const)
  call Control % Mass_Density     (dens_const)
  if(Flow % heat_transfer) then
    call Control % Heat_Capacity       (capa_const)
    call Control % Thermal_Conductivity(cond_const)
  end if
  if (Flow % n_scalars .gt. 0) then
    call Control % Scalars_Diffusivity (diff_const)
  end if

  if(Flow % with_interface) then
    call Control % Phase_Densities       (Vof % phase_dens)
    call Control % Phase_Viscosities     (Vof % phase_visc)
    call Control % Phase_Capacities      (Vof % phase_capa)
    call Control % Phase_Conductivities  (Vof % phase_cond)
    call Control % Surface_Tension       (Vof % surface_tension)
    call Control % Latent_Heat           (Vof % latent_heat)
    call Control % Saturation_Temperature(Vof % t_sat)
  else
    Flow % density       (:) = dens_const
    Flow % viscosity     (:) = visc_const
    if(Flow % heat_transfer) then
      Flow % capacity    (:) = capa_const
      Flow % conductivity(:) = cond_const
    end if
    if(Flow % n_scalars .gt. 0) then
      Flow % diffusivity (:) = diff_const
    end if
  end if

  end subroutine
