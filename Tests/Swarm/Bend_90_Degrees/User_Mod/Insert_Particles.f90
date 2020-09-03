!==============================================================================!
  subroutine User_Mod_Insert_Particles(flow, turb, mult, swarm, n, time)
!------------------------------------------------------------------------------!
!   This function is called at the end of time step.                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target     :: flow
  type(Turb_Type),       target     :: turb
  type(Multiphase_Type), target     :: mult
  type(Swarm_Type),      target     :: swarm
  integer,               intent(in) :: n         ! time step
  real,                  intent(in) :: time      ! physical time
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: grid
  integer :: k, n_parts_in_buffers, j, l, n_t, n_r, c, n_psr 
  real    :: dx, theta, theta_inc, sign_x, sign_z, rx, ry, rz
!==============================================================================!

  ! Take alias(es)
  grid => flow % pnt_grid

  !-------------------!
  !   1st time step   !
  !-------------------!

  n_psr     = int(sqrt(1. * swarm % n_particles))  ! sqrt(n_particles)
  n_r       = n_psr  ! number of nodes in azimuthal direction  (rotational)
  n_t       = n_psr  ! number of nodes in the radial direction (translational)
  theta_inc = 2.0 * PI / n_r

  if(n .eq. 2001) then     ! should be after the flow is developed

    ! Place the particles where you want them
    ! Theta loop
    do j = 1, n_r

      do k = 1, n_t

      ! particle index (1, 11, 21, 31, .... 2, 12, 22, 32, ....)
      l = (k - 1) * n_r + j

      ! Increasing theta (to swipe the area of the inlet)
      theta =  (j - 1) * theta_inc

      if(theta .ge. 0.0 .and. theta .lt. 90.0) then
        sign_x =  1.0
        sign_z = -1.0
      else if(theta .ge. 90.0 .and. theta .lt. 180.0) then 
        sign_x = -1.0
        sign_z = -1.0
      else if(theta .ge. 180.0 .and. theta .lt. 270.0) then
        sign_x = -1.0
        sign_z =  1.0
      else
        sign_x =  1.0
        sign_z =  1.0
      end if 

      ! Placing particles (only at the 1st time step)
      swarm % particle(l) % x_n = (0.009/n_psr) * cos(theta) * sign_x * k
      swarm % particle(l) % y_n =  0.138
      swarm % particle(l) % z_n = (0.009/n_psr) * sin(theta) * sign_z * k

      ! you essentially moved them a lot (from 0, 0, 0)
      swarm % particle(l) % cell = 0 
      swarm % particle(l) % node = 0 
      swarm % particle(l) % proc = 0 
      swarm % particle(l) % buff = 0

      swarm % particle(l) % x_o = swarm % particle(l) % x_n
      swarm % particle(l) % y_o = swarm % particle(l) % y_n
      swarm % particle(l) % z_o = swarm % particle(l) % z_n

      ! Searching for the closest cell and node to place the moved particle
      call Swarm_Mod_Find_Nearest_Cell(swarm, l, n_parts_in_buffers)
      call Swarm_Mod_Find_Nearest_Node(swarm, l)

      c = swarm % particle(l) % cell

      ! Set initial particle velocities
      rx = swarm % particle(l) % x_n - grid % xc(c)
      ry = swarm % particle(l) % y_n - grid % yc(c)
      rz = swarm % particle(l) % z_n - grid % zc(c)

      ! Compute velocities at the particle position from velocity gradients
      swarm % particle(l) % u    &
         = flow % u % n(c)       &  ! u velocity at the new time step (% n)
         + flow % u % x(c) * rx  &  ! u % x is gradient du/dx
         + flow % u % y(c) * ry  &  ! u % y is gradient du/dy
         + flow % u % z(c) * rz     ! u % x is gradient du/dz

      swarm % particle(l) % v    &
         = flow % v % n(c)       &  ! v velocity at the new time step (% n)
         + flow % v % x(c) * rx  &  ! v % x is gradient dv/dx
         + flow % v % y(c) * ry  &  ! v % y is gradient dv/dy
         + flow % v % z(c) * rz     ! v % x is gradient dv/dz

      swarm % particle(l) % w    &
         = flow % w % n(c)       &  ! w velocity at the new time step (% n)
         + flow % w % x(c) * rx  &  ! w % x is gradient dw/dx
         + flow % w % y(c) * ry  &  ! w % y is gradient dw/dy
         + flow % w % z(c) * rz     ! w % x is gradient dw/dz

      end do
    end do

  end if

  end subroutine
