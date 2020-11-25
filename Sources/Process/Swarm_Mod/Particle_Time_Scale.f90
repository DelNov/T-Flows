!==============================================================================!
  subroutine Swarm_Mod_Particle_Time_Scale(swarm)
!------------------------------------------------------------------------------!
!   Corrects particle time step size if input value is too big                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
  type(Turb_Type),  pointer :: turb
  real                      :: part_tss
  real                      :: visc_const
  real                      :: t_scale_max  ! max timescale for k-e-z-f model
  real                      :: t_ratio1
  real                      :: t_ratio2
  real                      :: t_ratio3
!==============================================================================!

  ! Take aliases for the swarm
  grid => swarm % pnt_grid
  flow => swarm % pnt_flow
  turb => swarm % pnt_turb

  ! Characteristic viscosity
  visc_const = maxval(flow % viscosity(:))

  ! Particle relaxation time
  swarm % tau = swarm % density * (swarm % diameter **2) / 18.0 / visc_const

  ! Calculate particle time step (tau_p/4)
  t_scale_max = maxval(turb % t_scale(:))
  part_tss    = swarm % tau * 0.25  ! probable particle time step size (Nyquist)
  t_ratio1    = flow % dt / part_tss     ! fluid tscale / part tscale critireon·
  t_ratio2    = t_scale_max / part_tss   ! turb tscale / part tscale critireon
  t_ratio3    = t_scale_max / flow % dt  ! turb tscale / flow tscale crit

  ! Correcting n_sup_steps (if the value in control file is small)
  if(t_ratio2 .ge. 1.0) then
    if(t_ratio1 .le. 1.0) then
      swarm % dt = flow % dt / real(swarm % n_sub_steps)
    else
      swarm % n_sub_steps = ceiling(t_ratio1)
      swarm % dt = flow % dt / real(swarm % n_sub_steps)
    end if
  else
    if(t_ratio3 .le. 1.0) then
      swarm % n_sub_steps = ceiling(t_ratio3)
      swarm % dt = t_scale_max
    else
      swarm % dt = flow % dt / real(swarm % n_sub_steps)
    end if
  end if

  end subroutine
