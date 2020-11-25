!==============================================================================!
  subroutine Swarm_Mod_Move_Particle(swarm, k)
!------------------------------------------------------------------------------!
!   Updates particle velocity and position                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
  integer                  :: k      ! particle number
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),    pointer :: flow
  type(Grid_Type),     pointer :: grid
  type(Turb_Type),     pointer :: turb
  type(Var_Type),      pointer :: u, v, w, kin, eps, zeta, f22
  type(Particle_Type), pointer :: part
  real,                pointer :: fb_x, fb_y, fb_z
  integer                      :: c, c2, n          ! nearest cell and node
  real                         :: dsp
  real                         :: rx, ry, rz        ! paticle-cell vector
  real                         :: dx, dy, dz        ! cell dimensions estimate
  real                         :: up, vp, wp        ! velocity at particle
  real                         :: flow_vel          ! flow vel. magnitude
  real                         :: k1, k2, k3, k4    ! Runge-Kutta increments
  real                         :: part_vel          ! relative velocity 
  real                         :: gravity           ! gravity magnitude
  real                         :: visc_const        ! characteristic viscosity
  real                         :: dens_const        ! characteristic density
  real                         :: f_fx, f_fy, f_fz  ! Brownian force components
  real                         :: fd_p              ! particle damping funct.
  real                         :: v2_mod_xc, v2_mod_yc, v2_mod_zc
  real                         :: u_mod, v_mod, w_mod
!==============================================================================!

  ! Take aliases for flow
  flow => swarm % pnt_flow
  grid => swarm % pnt_grid
  turb => swarm % pnt_turb
  call Field_Mod_Alias_Momentum(flow, u, v, w)

  ! Take aliases for turb
  call Turb_Mod_Alias_K_Eps_Zeta_F(turb, kin, eps, zeta, f22)

  ! Take aliases for swarm
  part => swarm % particle(k)
  fb_x => part  % fb_x
  fb_y => part  % fb_y
  fb_z => part  % fb_z

  ! Characteristic viscosity (needs to be discussed yet)
  ! (You can read it from the control file:
  ! call Control_Mod_Dynamic_Viscosity   (visc_const)
  ! call Control_Mod_Mass_Density        (dens_const)
  ! The way it is implemented now it could be different in every processor)
  visc_const = maxval(flow % viscosity(:))
  dens_const = maxval(flow % density(:))

  n  = part % node      ! index of the closest node for cfl calculation
  c  = part % cell      ! index of the closest cell for interpolation
  c2 = part % bnd_cell  ! index of the closest boundary cell for reflection

  ! Vector which connects particle position and cell centre
  rx = part % x_n - grid % xc(c)
  ry = part % y_n - grid % yc(c)
  rz = part % z_n - grid % zc(c)

  ! Particle damping function for wall treatment (for feeded modeled flow
  ! quantities)
  !fd_p = 1.0 - exp(-(turb % y_plus(c)/25.0)**3)     ! Piomelli
  if(turb % model .ne. NO_TURBULENCE_MODEL) then
    fd_p = (1.0 - exp(-(turb % y_plus(c)/25.0)))**2    ! Van-Driest
  end if

  !----------------------------------------------------------------------------!
  !   Compute velocities at the particle position from velocity gradients ...  !
  !   ... for all turbulence models, and even if there is no model at all      !
  !----------------------------------------------------------------------------!
  up = u % n(c)       &  ! u velocity at the new time step (% n)
     + u % x(c) * rx  &  ! u % x is gradient du/dx
     + u % y(c) * ry  &  ! u % y is gradient du/dy
     + u % z(c) * rz     ! u % x is gradient du/dz

  vp = v % n(c)       &  ! v velocity at the new time step (% n)
     + v % x(c) * rx  &  ! v % x is gradient dv/dx
     + v % y(c) * ry  &  ! v % y is gradient dv/dy
     + v % z(c) * rz     ! v % x is gradient dv/dz

  wp = w % n(c)       &  ! w velocity at the new time step (% n)
     + w % x(c) * rx  &  ! w % x is gradient dw/dx
     + w % y(c) * ry  &  ! w % y is gradient dw/dy
     + w % z(c) * rz     ! w % x is gradient dw/dz

  !--------------------------------------------!
  !   Add fluctuations from turbulent models   !
  !--------------------------------------------!
  if(turb % model .eq. HYBRID_LES_RANS) then

    ! Modeled quantity "zeta" at particle location
    v2_mod_xc = swarm % v2_mod_x(c) * rx
    v2_mod_yc = swarm % v2_mod_y(c) * ry
    v2_mod_zc = swarm % v2_mod_z(c) * rz

    ! Modeled wall-normal velocity component at particle location
    w_mod = swarm % v2_mod(c)  &
          + abs(v2_mod_xc)     &
          + abs(v2_mod_yc)     &
          + abs(v2_mod_zc)

!    ! Simpler version of the modeled quantities'
!    !   ad-hoc (No interpolation) "safer side"
!    if(turb % y_plus(c) > 35.0) then
!      w_mod = swarm % v2_mod(c) * fd_p
!    else 
!      w_mod = swarm % v2_mod(c)
!    end if 

    if(swarm % subgrid_scale_model .eq. DISCRETE_RANDOM_WALK) then

      ! Relative velocity components between particle and fluid
      part % rel_u_mod = abs(part % u - turb % u_mean(c))
      part % rel_v_mod = abs(part % v - turb % v_mean(c))
      part % rel_w_mod = abs(part % w - turb % w_mean(c))

      ! Compute the magnitude of the interpolated velocity
      flow_vel = sqrt(up**2 + vp**2 + wp**2)

      ! Compute the magnitude of the particle's velocity
      part_vel = sqrt(part % u **2 + part % v **2 + part % w **2)

      ! Compute Reynolds number for calculating Cd
      part % re = dens_const * part % d * abs(flow_vel - part_vel) / visc_const

      ! Compute the drag factor f
      part % f = 1.0 + 0.15 *(part % re ** 0.687)

      ! Add stochasticity by SEIM model (for modeled tubulent quantities)
      call Swarm_Mod_Sgs_Discrete_Random_Walk(swarm, k, rx, ry, rz)

      ! Adding the stochastic part from the Random Walk model (SEIM)
      up = up + part % u_drw
      vp = vp + part % v_drw
      wp = wp + part % v_drw

    else ! normal ER-HRL with ad-hoc SGS modeling 
      ! Compute wall-normal vel. at particle position from velocity gradients
      ! should create a switch for this case (SGS model is k-eps-zeta-f model
      ! itself)
      if(wp < 0.0) then
        w_mod = -1.0 * w_mod
      end if
      wp = wp + w_mod
    end if ! end of ER-HRL model with ad-hoc SGS modeling
  end if ! end of ER-HRL model with its two paths

  ! Compute the magnitude of the interpolated velocity
  flow_vel = sqrt(up**2 + vp**2 + wp**2)

  ! Compute the magnitude of the particle's velocity
  part_vel = sqrt(part % u **2 + part % v **2 + part % w **2)

  ! Particle relaxation time
  part % tau = part % density * (part % d **2) / 18.0 / visc_const

  ! Particle terminal speed 
  gravity = sqrt(grav_x**2 + grav_y**2 + grav_z**2)
  part % vel_t = part % tau * gravity

  ! Compute particle relative vel. in y-dir for buoyant force calculation
  part % rel_u   = up - part % u
  part % rel_v   = vp - part % v
  part % rel_w   = wp - part % w
  part % rel_vel = sqrt(  part % rel_u ** 2  &
                        + part % rel_v ** 2  &
                        + part % rel_w ** 2)

  ! Compute Reynolds number for calculating Cd
  part % re = dens_const * part % d * abs(flow_vel - part_vel) / visc_const

  ! Compute the drag factor f
  part % f = 1.0 + 0.15 *(part % re ** 0.687)

  !------------------------------------------!
  !                                          !
  !   Compute the new velocity of particle   !
  !           Runge-Kutta 4th order          !
  !                                          !
  !------------------------------------------!

  !----------------------------!
  !   Compute buoyancy force   !
  !----------------------------!
  fb_x = (1.0 - dens_const / part % density) * grav_x
  fb_y = (1.0 - dens_const / part % density) * grav_y
  fb_z = (1.0 - dens_const / part % density) * grav_z

  !------------------------!
  !   Fukagata SGS model   !
  !------------------------!
  f_fx = swarm % f_fuka_x(c)
  f_fy = swarm % f_fuka_y(c)
  f_fz = swarm % f_fuka_z(c)

  ! In the following we're accounting for only drag force & buoyancy force
  ! Attention: sign should be given to the gravity vector in the control file.
  ! ... according to the orientation relative to the global principal axes!
  !-------------------------!
  !   Updating x-velocity   !
  !-------------------------!
  k1 = (part % f * (up -  part % u) / part % tau) + fb_x + f_fx
  k2 = (part % f * (up - (part % u + (k1*swarm % dt)*0.5)) / part % tau)
  k3 = (part % f * (up - (part % u + (k2*swarm % dt)*0.5)) / part % tau)
  k4 = (part % f * (up - (part % u +  k3*swarm % dt))      / part % tau)

  ! X-velocity calculation
  part % u = part % u + (ONE_SIXTH) * (k1 + 2.0*(k2+k3) + k4)*swarm % dt

  !-------------------------!
  !   Updating y-velocity   !
  !-------------------------!
  k1 = (part % f * (vp -  part % v) / part % tau) + fb_y + f_fy
  k2 = (part % f * (vp - (part % v + (k1*swarm % dt)*0.5)) / part % tau)
  k3 = (part % f * (vp - (part % v + (k2*swarm % dt)*0.5)) / part % tau)
  k4 = (part % f * (vp - (part % v +  k3*swarm % dt))      / part % tau)

  ! Y-velocity calculation
  part % v = part % v + (ONE_SIXTH) * (k1 + 2.0*(k2+k3) + k4)*swarm % dt

  !-------------------------!
  !   Updating z-velocity   !
  !-------------------------!
  k1 = (part % f * (wp -  part % w) / part % tau) + fb_z + f_fz
  k2 = (part % f * (wp - (part % w + (k1*swarm % dt)*0.5)) / part % tau)
  k3 = (part % f * (wp - (part % w + (k2*swarm % dt)*0.5)) / part % tau)
  k4 = (part % f * (wp - (part % w +  k3*swarm % dt))      / part % tau)

  ! Z-velocity calculation
  part % w = part % w + (ONE_SIXTH) * (k1 + 2.0*(k2+k3) + k4)*swarm % dt

  !------------------------------------------------------------------!
  !   Compute the new position of particle with 1st order explicit   !
  !------------------------------------------------------------------!

  ! Update the particle position after reflection
  part % x_n = part % x_n + part % u * swarm % dt
  part % y_n = part % y_n + part % v * swarm % dt
  part % z_n = part % z_n + part % w * swarm % dt

  ! Calculating particle displacement
  dsp = Math_Mod_Distance_Squared(part % x_n, part % y_n, part % z_n, &
                                  part % x_o, part % y_o, part % z_o)

  ! Calculate cfl number for the particle (this is kind of approximate)
  dx = abs(grid % xn(n) - grid % xc(c)) * 2.0
  dy = abs(grid % yn(n) - grid % yc(c)) * 2.0
  dz = abs(grid % zn(n) - grid % zc(c)) * 2.0

  part % cfl = max(abs(part % u * swarm % dt) / dx,  &
                   abs(part % v * swarm % dt) / dy,  &
                   abs(part % w * swarm % dt) / dz)

  ! Particle stokes number
  ! St = tau_p/tau_f || tau_p = (rho_P*d_p^2)/18 mu || tau_f = nu/u_tau^2
  ! Should be done in a generic way by calling the friction velocity here...
  ! ... the used value used here is case_specific (Re_tau=590).
  part % st = swarm % density * (0.017046**2)       &
            * (part % d)**2 / 18.0 / visc_const**2

  end subroutine
