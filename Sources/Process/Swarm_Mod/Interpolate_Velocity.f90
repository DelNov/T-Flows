!==============================================================================!
  subroutine Swarm_Mod_Interpolate_Velocity(swarm, k)
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
  integer                      :: c                    ! nearest cell
  real                         :: rx, ry, rz           ! paticle-cell vector
  real                         :: x_old, y_old, z_old  ! particle's old x,y,z
  real                         :: dx, dy, dz           ! particle's shift
  real                         :: up, vp, wp           ! velocity at particle
  real                         :: flow_vel             ! flow vel. magn.
  real                         :: k1, k2, k3, k4       ! for Runge-Kutta
  real                         :: part_tau, part_cfl, part_vel
  logical,             pointer :: deposited            ! part. deposition flag
  logical,             pointer :: escaped              ! part. departure  flag
!==============================================================================!

  ! Take aliases
  flow      => swarm % pnt_flow
  grid      => swarm % pnt_grid
  u         => flow % u
  v         => flow % v
  w         => flow % w
  part      => swarm % particle(k)
  deposited => part  % deposited
  escaped   => part  % escaped

  c = part % cell ! assigning the index of the closest cell for interpolation

  ! Vector which connects particle position and cell centre
  rx = part % x - grid % xc(c)
  ry = part % y - grid % yc(c)
  rz = part % z - grid % zc(c)

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
  part % rel_vv = vp - part % v

  ! Compute Reynolds number for calculating Cd
  part % re = part % density * part % d * abs(flow_vel - part_vel) / viscosity

  ! Compute the drag factor f
  part % f = 1.0 + 0.15 *(part % re ** 0.687)

  !------------------------------------------!
  !   Compute the new velocity of particle   !
  !           Runge-Kutta 4th order          !
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

  !----------------------------------------!
  !  Compute the new position of particle  !
  !            1st order Explicit          !
  ! This step depends on the wall BC type  !
  !----------------------------------------!

  ! storing the old coordinates of particle before getting updated (for cfl)
  x_old = part % x
  y_old = part % y
  z_old = part % z

  !-----------------------------------!
  !    Trap condition (deposition)    !
  !-----------------------------------!
  if(swarm % rst <= TINY .and. .not. deposited) then

    ! Update the particle position after reflection
    part % x = part % x + part % u * swarm % dt
    part % y = part % y + part % v * swarm % dt
    part % z = part % z + part % w * swarm % dt

    ! Calculate cfl number for the particle  (trapped BC for walls)
    dx = abs(part % x - x_old)
    dy = abs(part % y - y_old)
    dz = abs(part % z - z_old)
    part_cfl = part_vel * swarm % dt / sqrt(dx**2 + dy**2 + dz**2)

    ! Printing particle position
    print *,k,'position','(',part%x,  &
    ',',part%y,',',part%z,')',',',' | cfl =',part_cfl

    if(part % y .le. 0.0) then    !just for the moment
      deposited = .true.
      swarm % cnt_d = swarm % cnt_d + 1
      print *, k, 'Particle is deposited!'
    end if
  end if

  !--------------------------!
  !   Reflection condition   !
  !--------------------------!
  if(swarm % rst > TINY) then
    if(part % y .le. 0.0000000000) then    !just for the moment until i make it generic
      part % y = 0.000001
      part % u = part % u * ( swarm % rst)
      part % v = part % v * (-swarm % rst)
      part % w = part % w * ( swarm % rst)

      ! Update the particle position after reflection
      part % x = part % x + part % u * swarm % dt
      part % y = part % y + part % v * swarm % dt
      part % z = part % z + part % w * swarm % dt

      ! Calculate cfl number for the particle (reflection BC)
      dx = abs(part % x - x_old)
      dy = abs(part % y - y_old)
      dz = abs(part % z - z_old)
      part_cfl = part_vel * swarm % dt / sqrt(dx**2 + dy**2 + dz**2)

      ! Increasing the number of particle reflections
      swarm % cnt_r = swarm % cnt_r + 1   ! to be engineered because ...
                                          ! ... a single particle can ...
                                          ! ... bounce several times.
      print *,k,'Particle is reflected!'
      print *,k,'position','(',part % x,  &
      ',',part % y,',',part % z,')',',',' | cfl =',part_cfl

      else
        ! If the particle didn't hit the wall, ...
        ! ... just update the position in the normal way
        part % x = part % x + part % u * swarm % dt
        part % y = part % y + part % v * swarm % dt
        part % z = part % z + part % w * swarm % dt

        ! Calculate cfl number for the particle (particles escape normally)
        dx = abs(part % x - x_old)
        dy = abs(part % y - y_old)
        dz = abs(part % z - z_old)
        part_cfl = (part_vel * swarm % dt) / sqrt(dx**2 + dy**2 + dz**2)

        print *,k,'position','(',part%x,  &
        ',',part%y,',',part%z,')',',',' | cfl =',part_cfl

      end if
   end if

  !--------------------------------------------------------------!
  !   Departure condition (particles escaping from the domain)   !
  !--------------------------------------------------------------!
  if(part % x .ge. 0.104999999999    &
  .or. part % x .le.  -0.104999999999) then    !just for the moment
    escaped =  .true.
    swarm % cnt_e = swarm % cnt_e + 1
    print *,k,'Particle escaped from outlet!'
  end if

  end subroutine
