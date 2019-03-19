!==============================================================================!
  subroutine Swarm_Mod_Move_Particle(swarm, k)
!------------------------------------------------------------------------------!
!   Interpolates velocity at the particle's position                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
  integer                  :: k      ! particle number
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),    pointer :: flow
  type(Grid_Type),     pointer :: grid
  type(Var_Type),      pointer :: u, v, w
  type(Particle_Type), pointer :: part
  integer                      :: c, c2                ! nearest cell
  real                         :: rx, ry, rz           ! paticle-cell vector
  real                         :: up, vp, wp           ! velocity at particle
  real                         :: flow_vel             ! flow vel. magn.
  real                         :: k1, k2, k3, k4       ! for Runge-Kutta
  real                         :: part_tau, part_vel
!==============================================================================!

  ! Take aliases
  flow => swarm % pnt_flow
  grid => swarm % pnt_grid
  u    => flow % u
  v    => flow % v
  w    => flow % w
  part => swarm % particle(k)

  c  = part % cell      ! index of the closest cell for interpolation
  c2 = part % bnd_cell  ! index of the closest boundary cell for reflection

  ! Vector which connects particle position and cell centre
  rx = part % x_n - grid % xc(c)
  ry = part % y_n - grid % yc(c)
  rz = part % z_n - grid % zc(c)

  ! Compute velocities at the particle position from velocity gradients
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

  ! Compute the magnitude of the interpolated velocity 
  flow_vel = sqrt(up**2 + vp**2 + wp**2)

  ! Compute the magnitude of the particle's velocity
  part_vel = sqrt(part % u **2 + part % v **2 + part % w **2)

  ! Particle relaxation time
  part_tau = part % density * (part % d **2) / 18.0 / viscosity

  ! Particle time step (division of the global time step)
  swarm % dt = flow % dt / 20.0

  ! Compute particle relative vel. in y-dir for buoyant force calculation
  part % rel_u   = up - part % u
  part % rel_v   = vp - part % v
  part % rel_w   = wp - part % w
  part % rel_vel = sqrt(  part % rel_u ** 2  &
                        + part % rel_v ** 2  &
                        + part % rel_w ** 2)

  ! Compute Reynolds number for calculating Cd
  part % re = part % density * part % d * abs(flow_vel - part_vel) / viscosity

  ! Compute the drag factor f
  part % f = 1.0 + 0.15 *(part % re ** 0.687)

  !------------------------------------------!
  !                                          !
  !   Compute the new velocity of particle   !
  !           Runge-Kutta 4th order          !
  !                                          !
  !------------------------------------------!

  !-------------------------!
  !   Updating x-velocity   !
  !-------------------------!
  k1 = part % f * (up -  part % u)                        / part_tau
  k2 = part % f * (up - (part % u + (k1*swarm % dt)*0.5)) / part_tau
  k3 = part % f * (up - (part % u + (k2*swarm % dt)*0.5)) / part_tau
  k4 = part % f * (up - (part % u +  k3*swarm % dt))      / part_tau

  ! X-velocity calculation
  part % u = part % u + (ONE_SIXTH) * (k1 + 2.0*(k2+k3) + k4)*swarm % dt

  !-------------------------!
  !   Updating y-velocity   !
  !-------------------------!
  k1 = (part % f * (vp -  part % v)                        / part_tau) - EARTH_G
  k2 = (part % f * (vp - (part % v + (k1*swarm % dt)*0.5)) / part_tau) - EARTH_G
  k3 = (part % f * (vp - (part % v + (k2*swarm % dt)*0.5)) / part_tau) - EARTH_G
  k4 = (part % f * (vp - (part % v +  k3*swarm % dt))      / part_tau) - EARTH_G

  ! Y-velocity calculation
  part % v = part % v + (ONE_SIXTH) * (k1 + 2.0*(k2+k3) + k4)*swarm % dt

  !-------------------------!
  !   Updating z-velocity   !
  !-------------------------!
  k1 = part % f * (wp -   part % w)                        / part_tau
  k2 = part % f * (wp -  (part % w + (k1*swarm % dt)*0.5)) / part_tau
  k3 = part % f * (wp -  (part % w + (k2*swarm % dt)*0.5)) / part_tau
  k4 = part % f * (wp -  (part % w +  k3*swarm % dt))      / part_tau

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

  ! Calculate cfl number for the particle (this is kind of wrong)
  part % cfl = part_vel * swarm % dt / grid % delta(c)

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
