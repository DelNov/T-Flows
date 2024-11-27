!==============================================================================!
  subroutine Particle_Time_Scale(Swarm)
!------------------------------------------------------------------------------!
!>  Adjusts the particle time step size within a swarm to ensure proper time
!>  scale resolution. This subroutine is pivotal for maintaining numerical
!>  stability and accuracy in simulations involving particle tracking, by
!>  adapting the time step size based on fluid and turbulence time scales.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Particle relaxation time: Computes the particle's characteristic time    !
!     scale based on its diameter and the fluid's viscosity.                   !
!   * Time scale ratios: Assesses the ratio of fluid and turbulence time       !
!     scales to the particle's time scale for guiding the adjustment of the    !
!     time step.                                                               !
!   * Time step adjustment: Modifies the number of sub-steps and particle      !
!     time step size to align with the Nyquist criterion and ensure adequate   !
!     resolution of particle dynamics relative to the flow field.              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Swarm_Type), target :: Swarm  !! the swarm of particles
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  type(Turb_Type),  pointer :: Turb
  real                      :: part_tss
  real                      :: visc_const
  real                      :: t_scale_max  ! max timescale for k-e-z-f model
  real                      :: t_ratio1
  real                      :: t_ratio2
  real                      :: t_ratio3
!==============================================================================!

  ! Take aliases for the swarm
  Grid => Swarm % pnt_grid
  Flow => Swarm % pnt_flow
  Turb => Swarm % pnt_turb

  ! Characteristic viscosity
  visc_const = maxval(Flow % viscosity(:))

  ! Particle relaxation time
  Swarm % tau = Swarm % density * (Swarm % diameter **2) / 18.0 / visc_const

  ! Calculate particle time step (tau_p/4)
  t_scale_max = maxval(Turb % t_scale(:))
  part_tss    = Swarm % tau * 0.25  ! probable particle time step size (Nyquist)
  t_ratio1    = Flow % dt / part_tss     ! fluid tscale / part tscale critireon·
  t_ratio2    = t_scale_max / part_tss   ! turb tscale / part tscale critireon
  t_ratio3    = t_scale_max / Flow % dt  ! turb tscale / Flow tscale crit

  ! Correcting n_sup_steps (if the value in control file is small)
  if(t_ratio2 .ge. 1.0) then
    if(t_ratio1 .le. 1.0) then
      Swarm % dt = Flow % dt / real(Swarm % n_sub_steps)
    else
      Swarm % n_sub_steps = ceiling(t_ratio1)
      Swarm % dt = Flow % dt / real(Swarm % n_sub_steps)
    end if
  else
    if(t_ratio3 .le. 1.0) then
      Swarm % n_sub_steps = ceiling(t_ratio3)
      Swarm % dt = t_scale_max
    else
      Swarm % dt = Flow % dt / real(Swarm % n_sub_steps)
    end if
  end if

  end subroutine
