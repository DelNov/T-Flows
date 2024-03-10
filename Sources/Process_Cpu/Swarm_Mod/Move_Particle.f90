!==============================================================================!
  subroutine Move_Particle(Swarm, k)
!------------------------------------------------------------------------------!
!>  Updates the position and velocity of a particle within a swarm. This
!>  subroutine is integral to the Lagrangian particle tracking mechanism in
!>  T-Flows. It computes new velocities and positions for each particle based
!>  on the current flow field, particle properties, and turbulence models.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Turbulence modeling: Integrates effects of turbulence on particle motion !
!     using subgrid scale (SGS) models like Brownian motion (Fukagata) and     !
!     Discrete Random Walk (DRW).                                              !
!   * Particle-flow interaction: Calculates particle velocities based on       !
!     local flow velocities, accounting for turbulence and fluid properties.   !
!   * Numerical integration: Employs a fourth-order Runge-Kutta method for     !
!     updating particle velocities and a first-order explicit scheme for       !
!     position updates.                                                        !
!   * Boundary interaction: Manages particle interaction with computational    !
!     domain boundaries, including reflection and wall effects.                !
!   * Physical forces: Considers forces like drag, buoyancy, and Brownian      !
!     motion in particle dynamics.                                             !
!   * Particle position Update: Determines new particle positions based on     !
!     updated velocities and time step size.                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Swarm_Type), target :: Swarm  !! the swarm of particles
  integer, intent(in)       :: k      !! particle number (rank)
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),    pointer :: Flow
  type(Grid_Type),     pointer :: Grid
  type(Turb_Type),     pointer :: Turb
  type(Var_Type),      pointer :: u, v, w, kin, eps, zeta, f22
  type(Particle_Type), pointer :: Part
  real,                pointer :: fb_x, fb_y, fb_z
  integer                      :: c, c2, n          ! nearest cell and node
  real                         :: dsp
  real                         :: rx, ry, rz        ! paticle-cell vector
  real                         :: dx, dy, dz        ! cell dimensions estimate
  real                         :: up, vp, wp        ! velocity at particle
  real                         :: flow_vel          ! Flow vel. magnitude
  real                         :: k1, k2, k3, k4    ! Runge-Kutta increments
  real                         :: part_vel          ! relative velocity 
  real                         :: gravity           ! gravity magnitude
  real                         :: visc_fluid        ! characteristic viscosity
  real                         :: dens_fluid        ! characteristic density
  real                         :: f_fx, f_fy, f_fz  ! Brownian force components
  real                         :: fd_p              ! particle damping funct.
  real                         :: v2_mod_xc, v2_mod_yc, v2_mod_zc
  real                         :: w_mod ! , u_mod, v_mod
!==============================================================================!

  ! Take aliases for Flow
  Flow => Swarm % pnt_flow
  Grid => Swarm % pnt_grid
  Turb => Swarm % pnt_turb
  call Flow % Alias_Momentum(u, v, w)

  ! Take aliases for Turb
  call Turb % Alias_K_Eps_Zeta_F(kin, eps, zeta, f22)

  ! Take aliases for Swarm
  Part => Swarm % Particle(k)
  fb_x => Part  % fb_x
  fb_y => Part  % fb_y
  fb_z => Part  % fb_z

  ! Take aliases for closest node, vell and boundary cell
  n  = Part % node      ! index of the closest node for cfl calculation
  c  = Part % cell      ! index of the closest cell for interpolation
  c2 = Part % bnd_cell  ! index of the closest boundary cell for reflection

  ! Characteristic viscosity and density (needs to be discussed yet)
  visc_fluid = Flow % viscosity(c)
  dens_fluid = Flow % density(c)

  ! Store fluid density for saving
  Part % dens_fluid = dens_fluid

  ! Vector which connects particle position and cell centre
  rx = Part % x_n - Grid % xc(c)
  ry = Part % y_n - Grid % yc(c)
  rz = Part % z_n - Grid % zc(c)

  ! Particle damping function for wall treatment (for feeded modeled Flow
  ! quantities)
  !fd_p = 1.0 - exp(-(Turb % y_plus(c)/25.0)**3)     ! Piomelli
  if(Turb % model .ne. NO_TURBULENCE_MODEL) then
    fd_p = (1.0 - exp(-(Turb % y_plus(c)/25.0)))**2    ! Van-Driest
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
  if(Turb % model .eq. HYBRID_LES_RANS) then

    ! Modeled quantity "zeta" at particle location
    v2_mod_xc = Swarm % v2_mod_x(c) * rx
    v2_mod_yc = Swarm % v2_mod_y(c) * ry
    v2_mod_zc = Swarm % v2_mod_z(c) * rz

    ! Modeled wall-normal velocity component at particle location
    w_mod = Swarm % v2_mod(c)  &
          + abs(v2_mod_xc)     &
          + abs(v2_mod_yc)     &
          + abs(v2_mod_zc)

!    ! Simpler version of the modeled quantities'
!    !   ad-hoc (No interpolation) "safer side"
!    if(Turb % y_plus(c) > 35.0) then
!      w_mod = Swarm % v2_mod(c) * fd_p
!    else 
!      w_mod = Swarm % v2_mod(c)
!    end if 

    if(Swarm % subgrid_scale_model .eq. DISCRETE_RANDOM_WALK) then

      ! Relative velocity components between particle and fluid
      Part % rel_u_mod = abs(Part % u - Turb % u_mean(c))
      Part % rel_v_mod = abs(Part % v - Turb % v_mean(c))
      Part % rel_w_mod = abs(Part % w - Turb % w_mean(c))

      ! Compute the magnitude of the interpolated velocity
      flow_vel = sqrt(up**2 + vp**2 + wp**2)

      ! Compute the magnitude of the particle's velocity
      part_vel = sqrt(Part % u **2 + Part % v **2 + Part % w **2)

      ! Compute Reynolds number for calculating Cd
      Part % re = dens_fluid * Part % d * abs(flow_vel - part_vel) / visc_fluid

      ! Compute the drag factor f
      Part % f = 1.0 + 0.15 *(Part % re ** 0.687)

      ! Add stochasticity by SEIM model (for modeled tubulent quantities)
      call Swarm % Sgs_Discrete_Random_Walk(k, rx, ry, rz)

      ! Adding the stochastic part from the Random Walk model (SEIM)
      up = up + Part % u_drw
      vp = vp + Part % v_drw
      wp = wp + Part % v_drw

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
  part_vel = sqrt(Part % u **2 + Part % v **2 + Part % w **2)

  ! Particle relaxation time
  Part % tau = Part % density * (Part % d **2) / 18.0 / visc_fluid

  ! Particle terminal speed 
  gravity = sqrt(Flow % grav_x**2 + Flow % grav_y**2 + Flow % grav_z**2)
  Part % vel_t = Part % tau * gravity

  ! Compute particle relative vel. in y-dir for buoyant force calculation
  Part % rel_u   = up - Part % u
  Part % rel_v   = vp - Part % v
  Part % rel_w   = wp - Part % w
  Part % rel_vel = sqrt(  Part % rel_u ** 2  &
                        + Part % rel_v ** 2  &
                        + Part % rel_w ** 2)

  ! Compute Reynolds number for calculating Cd
  Part % re = dens_fluid * Part % d * abs(flow_vel - part_vel) / visc_fluid

  ! Compute the drag factor f
  Part % f = 1.0 + 0.15 *(Part % re ** 0.687)

  !------------------------------------------!
  !                                          !
  !   Compute the new velocity of particle   !
  !           Runge-Kutta 4th order          !
  !                                          !
  !------------------------------------------!

  !----------------------------!
  !   Compute buoyancy force   !
  !----------------------------!
  fb_x = (1.0 - dens_fluid / Part % density) * Flow % grav_x
  fb_y = (1.0 - dens_fluid / Part % density) * Flow % grav_y
  fb_z = (1.0 - dens_fluid / Part % density) * Flow % grav_z

  !------------------------!
  !   Fukagata SGS model   !
  !------------------------!
  f_fx = Swarm % f_fuka_x(c)
  f_fy = Swarm % f_fuka_y(c)
  f_fz = Swarm % f_fuka_z(c)

  ! In the following we're accounting for only drag force & buoyancy force
  ! Attention: sign should be given to the gravity vector in the control file.
  ! ... according to the orientation relative to the global principal axes!
  !-------------------------!
  !   Updating x-velocity   !
  !-------------------------!
  k1 = (Part % f * (up -  Part % u) / Part % tau) + fb_x + f_fx
  k2 = (Part % f * (up - (Part % u + (k1*Swarm % dt)*0.5)) / Part % tau)
  k3 = (Part % f * (up - (Part % u + (k2*Swarm % dt)*0.5)) / Part % tau)
  k4 = (Part % f * (up - (Part % u +  k3*Swarm % dt))      / Part % tau)

  ! X-velocity calculation
  Part % u = Part % u + (ONE_SIXTH) * (k1 + 2.0*(k2+k3) + k4)*Swarm % dt

  !-------------------------!
  !   Updating y-velocity   !
  !-------------------------!
  k1 = (Part % f * (vp -  Part % v) / Part % tau) + fb_y + f_fy
  k2 = (Part % f * (vp - (Part % v + (k1*Swarm % dt)*0.5)) / Part % tau)
  k3 = (Part % f * (vp - (Part % v + (k2*Swarm % dt)*0.5)) / Part % tau)
  k4 = (Part % f * (vp - (Part % v +  k3*Swarm % dt))      / Part % tau)

  ! Y-velocity calculation
  Part % v = Part % v + (ONE_SIXTH) * (k1 + 2.0*(k2+k3) + k4)*Swarm % dt

  !-------------------------!
  !   Updating z-velocity   !
  !-------------------------!
  k1 = (Part % f * (wp -  Part % w) / Part % tau) + fb_z + f_fz
  k2 = (Part % f * (wp - (Part % w + (k1*Swarm % dt)*0.5)) / Part % tau)
  k3 = (Part % f * (wp - (Part % w + (k2*Swarm % dt)*0.5)) / Part % tau)
  k4 = (Part % f * (wp - (Part % w +  k3*Swarm % dt))      / Part % tau)

  ! Z-velocity calculation
  Part % w = Part % w + (ONE_SIXTH) * (k1 + 2.0*(k2+k3) + k4)*Swarm % dt

  !------------------------------------------------------------------!
  !   Compute the new position of particle with 1st order explicit   !
  !------------------------------------------------------------------!

  ! Update the particle position after reflection
  Part % x_n = Part % x_n + Part % u * Swarm % dt
  Part % y_n = Part % y_n + Part % v * Swarm % dt
  Part % z_n = Part % z_n + Part % w * Swarm % dt

  ! Calculating particle displacement
  dsp = Math % Distance_Squared(Part % x_n, Part % y_n, Part % z_n, &
                                Part % x_o, Part % y_o, Part % z_o)

  ! Calculate cfl number for the particle (this is kind of approximate)
  dx = abs(Grid % xn(n) - Grid % xc(c)) * 2.0
  dy = abs(Grid % yn(n) - Grid % yc(c)) * 2.0
  dz = abs(Grid % zn(n) - Grid % zc(c)) * 2.0

  Part % cfl = max(abs(Part % u * Swarm % dt) / dx,  &
                   abs(Part % v * Swarm % dt) / dy,  &
                   abs(Part % w * Swarm % dt) / dz)

  ! Particle stokes number
  ! St = tau_p/tau_f || tau_p = (rho_P*d_p^2)/18 mu || tau_f = nu/u_tau^2
  ! Should be done in a generic way by calling the friction velocity here...
  ! ... the used value used here is case_specific (Re_tau=590).
  Part % st = Swarm % density * (0.017046**2)       &
            * (Part % d)**2 / 18.0 / visc_fluid**2

  end subroutine
