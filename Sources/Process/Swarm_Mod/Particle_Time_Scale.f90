!==============================================================================!
  subroutine Swarm_Mod_Particle_Time_Scale(swarm, turb)
!------------------------------------------------------------------------------!
!         Corrects particle time step size if input value is too big           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
  type(Turb_Type),  target :: turb
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: grid
  type(Field_Type),    pointer :: flow
  real                         :: part_tss      ! particle timestep size criterion
  real                         :: flow_dt       ! flow global time step
  real                         :: visc_const    ! C/c viscosity 
  real                         :: t_scale_max   ! maximum timescale for k-e-z-f model
  real                         :: t_ratio1      ! fluid time scale / particle_tscale_critireon 
  real                         :: t_ratio2      ! turb_timescale / particle_tscale_critireon
  real                         :: t_ratio3      ! turb_timescale / flow global time step
!==============================================================================!

  ! Take aliases for the swarm
  grid => swarm % pnt_grid
  flow => swarm % pnt_flow

  ! Characteristic viscosity 
  visc_const = maxval(flow % viscosity(:))

  ! Particle relaxation time
  swarm % tau = swarm % density * (swarm % diameter **2) / 18.0 / visc_const

  ! Calculate particle time step (tau_p/4)
  t_scale_max = maxval(turb % t_scale(:))
  part_tss    = swarm % tau * 0.25  ! probable particle time step size (Nyquist)
  flow_dt     = flow % dt
  t_ratio1    = flow_dt / part_tss 
  t_ratio2    = t_scale_max / (part_tss) 
  t_ratio3    = t_scale_max / flow_dt 

  ! correcting n_sup_steps (if the value in control file is small)
  if(t_ratio2 .ge. 1.0) then
    if(t_ratio1 .le. 1.0) then 
      !swarm % n_sub_steps = 50 
      swarm % dt = flow_dt / swarm % n_sub_steps
    else 
      swarm % n_sub_steps = ceiling(t_ratio1)
      swarm % dt = flow_dt / swarm % n_sub_steps 
    end if  
  else 
    if(t_ratio3 .le. 1.0) then
      swarm % n_sub_steps = ceiling(t_ratio3)
      swarm % dt = t_scale_max
    else 
      !swarm % n_sub_steps = 50
      swarm % dt = flow_dt / swarm % n_sub_steps
    end if 
  end if 

  end subroutine
