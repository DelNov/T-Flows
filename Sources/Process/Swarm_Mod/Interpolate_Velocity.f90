!==============================================================================!
  subroutine Swarm_Mod_Interpolate_Velocity(flow, swarm, k)
!------------------------------------------------------------------------------!
!               Interpolates velocity at the particle's position               !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod,      only: Field_Type, viscosity
  use Grid_Mod,       only: Grid_Type
  use Var_Mod,        only: Var_Type
  use Control_Mod
  use Const_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: flow
  type(Swarm_Type), target    :: swarm
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: grid
  type(Var_Type),      pointer :: u, v, w
  type(Particle_Type), pointer :: part
  integer                  :: c, n                 ! nearest cell and node
  integer                  :: n_part               ! number of particles
  integer                  :: k                    ! counter for particles
  real                     :: rx, ry, rz           ! vector connecting paticle and cell
  real                     :: xc, yc, zc           ! nearest cell coordinates
  real                     :: xn, yn, zn           ! nearest node coordinates
  real                     :: x_old, y_old, z_old  ! particle's old coordinates
  real                     :: dx, dy, dz           ! particle length scale in 3D
  real                     :: dr                   ! particle distance vectorfor cfl calc.
  real                     :: up, vp, wp           ! velocity at particle position
  real                     :: u_f                  ! flow velocity magnitude 
  real                     :: vel_p                ! particle velocity magnitude
  real                     :: f                    ! drag factor for calculating Cd 
  real                     :: Re                   ! Reynolds number (Particle)
  real                     :: dt                   ! time step
  real                     :: k1, k2, k3, k4       ! Runge-Kutta coefficients
  real                     :: test1, test2, test3  ! for debugguing
  real                     :: st                   ! stokes number
  real                     :: rst                  ! Coeff. of restitution
  real                     :: diameter             ! pipe diameter
  logical,             pointer :: deposited            ! particle deposition flag
  logical,             pointer:: escaped              ! particle departure  flag
  logical,             pointer:: trapped              ! trap       BC type
  logical,             pointer:: reflected            ! reflection BC type
!==============================================================================!

  ! Take aliases
  grid         => flow % pnt_grid
  u            => flow % u
  v            => flow % v
  w            => flow % w
  part         => swarm % particles(k)
  deposited    => part % deposited
  escaped      => part % escaped
  trapped      => part % trapped
  reflected    => part % reflected

  ! Importing the global time step (this slows down a bit)
  call Control_Mod_Time_Step(dt, verbose=.true.)

  c = part % cell ! assigning the index of the closest cell for interpolation 

  ! Cell centre coorindates
  xc = grid % xc(c)
  yc = grid % yc(c)
  zc = grid % zc(c)

  ! Closest node coordinates
  xn = grid % xn(n)
  yn = grid % yn(n)
  zn = grid % zn(n)

  ! Vector which connects particle position and cell centre
  rx = part % x - xc
  ry = part % y - yc
  rz = part % z - zc

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
  u_f = sqrt(up**2 + vp**2 + wp**2)

  ! Compute the magnitude of the particle's velocity
  part %vel_p = sqrt(part % u **2 + part % v **2 + part % w **2 )

  ! Particle's characteristics
  part % rho = 1000.0
  part % d   = 2.5e-5        !particle diameter

  ! Coefficient of restitution (particle reflection)
  rst = 1.0

  ! Pipe diameter
  diameter = 0.01

  ! Particle relaxation time
  part%tau = part % rho * (part % d **2) / 18.0 / viscosity

  ! Compute stokes number
  St= u_f * part%tau/ diameter

  ! Particle time step (division of the global time step)
  !part % dtp  = dt / 100.0
  part % dtp  = dt / 20.0

  ! Compute particle relative vel. in y-dir for buoyant force calculation
  part % rel_vv   =  vp - part % v
  part % rel_vvn  =  abs(vp - part % v)

  ! Compute Reynolds number for calculating Cd
  part % rel_vn   = abs(u_f - part % vel_p)
  part % Re = part % rho * part % d * (part %rel_vn) / viscosity

  ! Compute the drag factor f
  part % f = 1.0 + 0.15 *(part % Re **0.687)


  !------------------------------------------!
  !   Compute the new velocity of particle   !
  !           Runge-Kutta 4th order          !
  !------------------------------------------!

  !-------------------------!
  !   Updating x-velocity   !
  !-------------------------!
  k1 = part % f * (up -  part % u)                      / part % tau
  k2 = part % f * (up - (part % u + (k1*part%dtp)*0.5)) / part % tau
  k3 = part % f * (up - (part % u + (k2*part%dtp)*0.5)) / part % tau
  k4 = part % f * (up - (part % u +  k3*part%dtp))      / part % tau

  ! X-velocity calculation
  part % u = part % u + (ONE_SIXTH) * (k1 + 2.0*(k2+k3) + k4)*part % dtp

  !-------------------------!
  !   Updating y-velocity   !
  !-------------------------!
  k1 = (part % f * (vp -  part % v)                      / part % tau) - EARTH_G
  k2 = (part % f * (vp - (part % v + (k1*part%dtp)*0.5)) / part % tau) - EARTH_G
  k3 = (part % f * (vp - (part % v + (k2*part%dtp)*0.5)) / part % tau) - EARTH_G
  k4 = (part % f * (vp - (part % v +  k3*part%dtp))      / part % tau) - EARTH_G

  ! Y-velocity calculation
  part % v = part % v + (ONE_SIXTH) * (k1 + 2.0*(k2+k3) + k4)*part%dtp

  !-------------------------!
  !   Updating z-velocity   !
  !-------------------------!
  k1 = part % f * (wp -   part % w)                      / part % tau
  k2 = part % f * (wp -  (part % w + (k1*part%dtp)*0.5)) / part % tau
  k3 = part % f * (wp -  (part % w + (k2*part%dtp)*0.5)) / part % tau
  k4 = part % f * (wp -  (part % w +  k3*part%dtp))      / part % tau

  ! Z-velocity calculation
  part % w = part % w + (ONE_SIXTH) * (k1 + 2.0*(k2+k3) + k4)*part%dtp

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

  if (trapped) then
    ! Update the particle position after reflection
    part % x = part % x + part % u * part%dtp
    part % y = part % y + part % v * part%dtp
    part % z = part % z + part % w * part%dtp

    ! Calculate cfl number for the particle  (trapped BC for walls)
    dx = abs(part % x - x_old)
    dy = abs(part % y - y_old)
    dz = abs(part % z - z_old)
    dr = sqrt(dx**2 + dy**2 + dz**2)
    part % cfl = part %vel_p * part % dtp /dr

    ! Printing particle position
    print *,k,'position','(',part%x,  &
    ',',part%y,',',part%z,')',',',' | cfl =',part % cfl

    if(part % y .le. 0.0) then    !just for the moment
      deposited =  .true.
      swarm% c_d = swarm% c_d + 1
      print *,k,'Particle is deposited!'
    end if
  end if

  !--------------------------!
  !   Reflection condition   !
  !--------------------------!
  if (reflected) then
    if(part % y .le. 0.0000000000) then    !just for the moment until i make it generic
      part%y=0.000001
      part%u=part%u * ( rst)
      part%v=part%v * (-rst)
      part%w=part%w * ( rst)

      ! Update the particle position after reflection
      part % x = part % x + part % u * part%dtp
      part % y = part % y + part % v * part%dtp
      part % z = part % z + part % w * part%dtp

      ! Calculate cfl number for the particle (reflection BC)
      dx = abs(part % x - x_old)
      dy = abs(part % y - y_old)
      dz = abs(part % z - z_old)
      dr = sqrt(dx**2 + dy**2 + dz**2)
      part % cfl = part %vel_p * part % dtp /dr

      ! Increasing the number of particle reflections
      swarm% c_r = swarm% c_r + 1          ! To be engineered because a single particle can bounce several times..
      print *,k,'Particle is reflected!'
      print *,k,'position','(',part%x,  &
      ',',part%y,',',part%z,')',',',' | cfl =',part % cfl

      else
        ! If the particle didn't hit the wall, just update the position in the normal way
        part % x = part % x + part % u * part%dtp
        part % y = part % y + part % v * part%dtp
        part % z = part % z + part % w * part%dtp

        ! Calculate cfl number for the particle (particles escape normally)
        dx = abs(part % x - x_old)
        dy = abs(part % y - y_old)
        dz = abs(part % z - z_old)
        dr = sqrt(dx**2 + dy**2 + dz**2)
        part % cfl = (part %vel_p * part % dtp) / dr

        print *,k,'position','(',part%x,  &
        ',',part%y,',',part%z,')',',',' | cfl =',part % cfl

      end if
   end if

  !--------------------------------------------------------------!
  !   Departure condition (particles escaping from the domain)   !
  !--------------------------------------------------------------!
  if(part % x .ge. 0.104999999999    &
  .or. part % x .le.  -0.104999999999) then    !just for the moment
    escaped =  .true.
    swarm%c_e = swarm% c_e + 1
    print *,k,'Particle escaped from outlet!'
  end if

  end subroutine



