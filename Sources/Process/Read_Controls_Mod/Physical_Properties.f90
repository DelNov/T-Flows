!==============================================================================!
  subroutine Physical_Properties(Rc, Flow, Vof, Swarm)
!------------------------------------------------------------------------------!
!   Reads physical properties from control file.                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Read_Controls_Type), intent(in) :: Rc
  type(Field_Type)                      :: Flow
  type(Vof_Type)                        :: Vof
  type(Swarm_Type)                      :: Swarm
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  real                     :: dens_const, visc_const
  real                     :: capa_const, cond_const
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Rc)
  Unused(Swarm)
!==============================================================================!

  ! Give some sign
  if(First_Proc())  &
    print '(a)', ' # Reading about physical properties'

  ! Take alias
  Grid => Flow % pnt_grid

  ! Read constant (defualt) values
  call Control % Dynamic_Viscosity   (visc_const)
  call Control % Mass_Density        (dens_const)
  call Control % Heat_Capacity       (capa_const)
  call Control % Thermal_Conductivity(cond_const)
  call Control % Scalars_Diffusivity (Flow % diffusivity)

  if(Flow % with_interface) then
    call Control % Phase_Densities       (Vof % phase_dens)
    call Control % Phase_Viscosities     (Vof % phase_visc)
    call Control % Phase_Capacities      (Vof % phase_capa)
    call Control % Phase_Conductivities  (Vof % phase_cond)
    call Control % Surface_Tension       (Vof % surface_tension)
    call Control % Latent_Heat           (Vof % latent_heat)
    call Control % Saturation_Temperature(Vof % t_sat)
  else
    Flow % density     (:) = dens_const
    Flow % viscosity   (:) = visc_const
    Flow % capacity    (:) = capa_const
    Flow % conductivity(:) = cond_const
  end if

  end subroutine
