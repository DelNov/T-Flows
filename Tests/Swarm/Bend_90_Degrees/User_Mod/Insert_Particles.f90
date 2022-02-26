!==============================================================================!
  subroutine User_Mod_Insert_Particles(Flow, Turb, Vof, Swarm, n, time)
!------------------------------------------------------------------------------!
!   This function is called at the end of time step.                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
  integer, intent(in)      :: n       ! time step
  real,    intent(in)      :: time    ! physical time
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: i, j, k, l, n_parts_in_buffers, n_r, c
  real                     :: rx, ry, rz, r, theta, d_r, d_theta, x, z
!==============================================================================!

  ! Take alias(es)
  Grid => Flow % pnt_grid

  n_r = 40

  !-------------------!
  !   1st time step   !
  !-------------------!
  if(n .eq. 2401) then     ! should be after the Flow is developed

    ! First Particle in the center
    l = 1
    Swarm % Particle(l) % x_n =  0.0
    Swarm % Particle(l) % y_n =  0.0399999
    Swarm % Particle(l) % z_n =  0.0

    d_r = 0.0085 / (n_r - 2)

    ! Place the particles where you want them
    ! Theta loop
    do i = 1, n_r - 1
      r = i * d_r
      d_theta = 2.0 * PI / (i * 6)
      do j = 1, i * 6
        theta = (j-1) * d_theta
        x = r * cos(theta)
        z = r * sin(theta)

        l = l + 1
        Swarm % Particle(l) % x_n = x
        Swarm % Particle(l) % y_n = 0.0399
        Swarm % Particle(l) % z_n = z
      end do
    end do

    if(l <= 10000) then
      if(this_proc < 2) then
        print *, '# @User_Mod_Insert_Particles: inserted', l, ' particles'
      end if
      Swarm % n_particles = l
    else
      if(this_proc < 2) then
        print *, '# @User_Mod_Insert_Particles: too many patrticles, reduce n_r'
        call Comm_Mod_End
        stop
      end if
    end if

    do l = 1, Swarm % n_particles

      ! You essentially moved them a lot (from 0, 0, 0)
      Swarm % Particle(l) % cell = 0 
      Swarm % Particle(l) % node = 0 
      Swarm % Particle(l) % proc = 0 
      Swarm % Particle(l) % buff = 0

      Swarm % Particle(l) % x_o = Swarm % Particle(l) % x_n
      Swarm % Particle(l) % y_o = Swarm % Particle(l) % y_n
      Swarm % Particle(l) % z_o = Swarm % Particle(l) % z_n

      ! Searching for the closest cell and node to place the moved Particle
      call Swarm % Particle(l) % Find_Nearest_Cell(n_parts_in_buffers)
      call Swarm % Particle(l) % Find_Nearest_Node()

      c = Swarm % Particle(l) % cell

      ! Set initial Particle velocities
      rx = Swarm % Particle(l) % x_n - Grid % xc(c)
      ry = Swarm % Particle(l) % y_n - Grid % yc(c)
      rz = Swarm % Particle(l) % z_n - Grid % zc(c)

      ! Compute velocities at the Particle position from velocity gradients
      Swarm % Particle(l) % u    &
         = Flow % u % n(c)       &  ! u velocity at the new time step (% n)
         + Flow % u % x(c) * rx  &  ! u % x is gradient du/dx
         + Flow % u % y(c) * ry  &  ! u % y is gradient du/dy
         + Flow % u % z(c) * rz     ! u % x is gradient du/dz

      Swarm % Particle(l) % v    &
         = Flow % v % n(c)       &  ! v velocity at the new time step (% n)
         + Flow % v % x(c) * rx  &  ! v % x is gradient dv/dx
         + Flow % v % y(c) * ry  &  ! v % y is gradient dv/dy
         + Flow % v % z(c) * rz     ! v % x is gradient dv/dz

      Swarm % Particle(l) % w    &
         = Flow % w % n(c)       &  ! w velocity at the new time step (% n)
         + Flow % w % x(c) * rx  &  ! w % x is gradient dw/dx
         + Flow % w % y(c) * ry  &  ! w % y is gradient dw/dy
         + Flow % w % z(c) * rz     ! w % x is gradient dw/dz

    end do

  end if

  end subroutine
