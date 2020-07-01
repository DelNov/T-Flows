!==============================================================================!
  subroutine Swarm_Mod_Move_Particle(swarm, turb, k)
!------------------------------------------------------------------------------!
!   Updates particle velocity and position                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
  type(Turb_Type),  target :: turb
  integer                  :: k      ! particle number
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),    pointer :: flow
  type(Grid_Type),     pointer :: grid
  type(Var_Type),      pointer :: u, v, w, kin, eps, zeta, f22
  type(Particle_Type), pointer :: part
  real,                pointer :: fb_x, fb_y, fb_z
  integer                      :: c, c2             ! nearest cell
  real                         :: rx, ry, rz        ! paticle-cell vector
  real                         :: up, vp, wp        ! velocity at particle
  real                         :: flow_vel          ! flow vel. magnitude
  real                         :: k1, k2, k3, k4    ! Runge-Kutta increments
  real                         :: part_vel          ! relative velocity 
  real                         :: gravity           ! gravity magnitude
  real                         :: visc_const        ! characteristic viscosity
  real                         :: dens_const        ! characteristic density
  real                         :: f_fx, f_fy, f_fz  ! Brownian force components
  real                         :: fd_p              ! particle damping fn. 
  real                         :: t_scale_max       ! maximum t_scale for ERHRL 
                                                    ! model on the grid used
  ! Gradients of modeled flow velocities
  real                         :: u_mod_xc, v_mod_xc, w_mod_xc
  real                         :: u_mod_yc, v_mod_yc, w_mod_yc
  real                         :: u_mod_zc, v_mod_zc, w_mod_zc
!==============================================================================!

  ! Take aliases for flow
  flow => swarm % pnt_flow
  grid => swarm % pnt_grid
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

  c  = part % cell      ! index of the closest cell for interpolation
  c2 = part % bnd_cell  ! index of the closest boundary cell for reflection

  ! Vector which connects particle position and cell centre
  rx = part % x_n - grid % xc(c)
  ry = part % y_n - grid % yc(c)
  rz = part % z_n - grid % zc(c)

  ! Particle damping function for wall treatment (for feeded modeled flow
  ! quantities)
  fd_p = 1.0 - exp(-(turb % y_plus(c)/25.0)**3)

  ! Compute velocities at the particle position from velocity gradients
  if(turb % model .eq. HYBRID_LES_PRANDTL) then
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
  end if

  if(turb % model .eq. HYBRID_LES_RANS) then

    ! Modeled quantity "zeta" at particle location  
    w_mod_xc = sqrt(swarm % w_mod_x(c) * rx * rx)  
    w_mod_yc = sqrt(swarm % w_mod_y(c) * ry * ry)   
    w_mod_zc = sqrt(swarm % w_mod_z(c) * rz * rz)

    ! The proportion between u' and v' is rudimentarily computed from DNS data
    ! at
    ! Re=180 in periodic channel flow (just as a rough estimation at the moment)
    ! ATTENTION! we don't have the gradients in u and v directions...!
    u_mod_xc = (kin % n(c) - w_mod_xc) / 1.2757
    u_mod_yc = (kin % n(c) - w_mod_yc) / 1.2757
    u_mod_zc = (kin % n(c) - w_mod_zc) / 1.2757

    v_mod_xc = 0.2757 * u_mod_xc
    v_mod_yc = 0.2757 * u_mod_yc
    v_mod_zc = 0.2757 * u_mod_zc

    ! Compute velocities at the particle position from velocity gradients
    ! should create a switch for this case (SGS model is k-eps-zeta-f model
    ! itself)
    up = u % n(c)       &  ! u velocity at the new time step (% n)
       + u % x(c) * rx  &  ! u % x is gradient du/dx
       + u % y(c) * ry  &  ! u % y is gradient du/dy
       + u % z(c) * rz  &  ! u % x is gradient du/dz
       + u_mod_xc       &  ! u_mod_dx is gradient du_mod/dx 
       + u_mod_yc       &  ! u_mod_dy is gradient du_mod/dy   
       + u_mod_zc          ! u_mod_dz is gradient du_mod/dz  

    vp = v % n(c)      &  ! v velocity at the new time step (% n)
       + v % x(c) * rx  &  ! v % x is gradient dv/dx
       + v % y(c) * ry  &  ! v % y is gradient dv/dy
       + v % z(c) * rz  &  ! v % x is gradient dv/dz
       + v_mod_xc       &  ! v_mod_dx is gradient dv_mod/dx
       + v_mod_yc       &  ! v_mod_dy is gradient dv_mod/dy
       + v_mod_zc          ! v_mod_dz is gradient dv_mod/dz

    wp = w % n(c)       &      ! w velocity at the new time step (% n)
       + w % x(c) * rx  &      ! w % x is gradient dw/dx
       + w % y(c) * ry  &      ! w % y is gradient dw/dy
       + w % z(c) * rz  &      ! w % x is gradient dw/dz
       + w_mod_xc       &      ! w_mod_dx is gradient dw_mod/dx
       + w_mod_yc       &      ! w_mod_dy is gradient dw_mod/dy
       + w_mod_zc              ! w_mod_dz is gradient dw_mod/dz
  end if

  ! Compute the magnitude of the interpolated velocity 
  flow_vel = sqrt(up**2 + vp**2 + wp**2)

  ! Compute the magnitude of the particle's velocity
  part_vel = sqrt(part % u **2 + part % v **2 + part % w **2)

  ! Particle relaxation time
  part % tau = part % density * (part % d **2) / 18.0 / visc_const

  ! Particle terminal speed 
  GRAVITY = sqrt(grav_x**2 + grav_y**2 + grav_z**2)
  part % vel_t = part % tau * GRAVITY 

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
  !if(this_proc .le. 2) then
  !  call Sgs_Fukagata_Model(swarm, turb, k)
  !end if 

  f_fx = swarm % f_fuka_x(c)
  f_fy = swarm % f_fuka_y(c)
  f_fz = swarm % f_fuka_z(c)

  ! In the following we're accounting for only drag force & buoyancy force
  ! Attention: sign should be given to the gravity vector in the control file...
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

  ! storing the old coordinates of particle before getting updated (for cfl)
  part % x_o = part % x_n
  part % y_o = part % y_n
  part % z_o = part % z_n

  ! Update the particle position after reflection
  part % x_n = part % x_n + part % u * swarm % dt
  part % y_n = part % y_n + part % v * swarm % dt
  part % z_n = part % z_n + part % w * swarm % dt

  ! Calculating the distance that particle just moved (for cfl condition)
  part % dsp = Math_Mod_Distance_Squared(part % x_n, part % y_n, part % z_n, &
                                         part % x_o, part % y_o, part % z_o)

  ! Calculate cfl number for the particle (this is kind of wrong)
  part % cfl = part_vel * swarm % dt / part % dsp

  ! Particle stokes number 
  ! St= tau_p/tau_f || tau_p = (rho_P*d_p^2)/18 mu || tau_f = nu/u_tau^2    ...
  !... should be done in a generic way by calling the friction velocity here...
  !... the used value used here is case_specific (Re_tau=150).
  swarm % st = swarm % density * (part % d)**2/18.0/(visc_const**2)*0.0043**2

  !----------------------------------------------!
  !                                              !
  !   If particle is close to a boundary cell,   !
  !    handle its interaction with boundaries    !
  !                                              !
  !----------------------------------------------!
  if(c2 .ne. 0) then
    call Swarm_Mod_Bounce_Particle(swarm, k)
  end if

  end subroutine
