!==============================================================================!
  subroutine User_Mod_Insert_Particles(Flow, Turb, Vof, Swarm, curr_dt, time)
!------------------------------------------------------------------------------!
!   This function is called at the end of time step.                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
  integer, intent(in)      :: curr_dt  ! time step
  real,    intent(in)      :: time     ! physical time
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: i, j, k, c, n_parts_in_buffers
  integer                  :: n_rows, n_new, n_old
  real                     :: rx, ry, rz, r, theta, d_r, d_theta, x, z
!==============================================================================!

  ! Take alias(es)
  Grid => Flow % pnt_grid

  ! Re-set the diameters, since what you have in backup could be different
  Swarm % Particle(:) % d =  Swarm % diameter

  !----------------------------------------!
  !   Number of rows in radial direction   !
  !----------------------------------------!

  ! Number of particles = 1 + 6 + 12 + ... (n_rows - 1) * 6
  !                     = 1 + 6 * (1 + 2 + 3 + ... + n_rows)
  !                     = 1 + 3 * n_rows * (n_rows-1)
  n_rows = 15

  !-------------------------------------------!
  !                                           !
  !   1st time step of particle computation   !
  !                                           !
  !-------------------------------------------!
  n_old = Swarm % n_particles   ! old number of particles

  if(curr_dt .eq. 1001 .or.  &
     curr_dt .eq. 1501 .or.  &
     curr_dt .eq. 2001 .or.  &
     curr_dt .eq. 2501) then      ! should be after the flow is developed

    ! First Particle in the center
    k = n_old + 1
    Swarm % Particle(k) % x_n =  0.0
    Swarm % Particle(k) % y_n =  0.0999
    Swarm % Particle(k) % z_n =  0.0

    d_r = 0.009 / (n_rows - 2)

    ! Place the particles where you want them
    ! Theta loop
    do i = 1, n_rows - 1
      r = i * d_r
      d_theta = 2.0 * PI / (i * 6)
      do j = 1, i * 6
        theta = (j-1) * d_theta
        x = r * cos(theta)
        z = r * sin(theta)

        k = k + 1
        Swarm % Particle(k) % x_n = x
        Swarm % Particle(k) % y_n = 0.0999
        Swarm % Particle(k) % z_n = z
      end do
    end do

    n_new = k  ! new number of particles

    if(n_new <= Swarm % max_particles) then
      if(this_proc < 2) then
        print *, '# @User_Mod_Insert_Particles: inserted',  &
                 n_new - n_old, ' particles'
      end if
    else
      if(this_proc < 2) then
        print *, '# @User_Mod_Insert_Particles: too many particles'
        call Comm_Mod_End
        stop
      end if
    end if

    do k = n_old + 1, n_new

      ! You essentially moved them a lot (from 0, 0, 0)
      Swarm % Particle(k) % cell = 0
      Swarm % Particle(k) % node = 0
      Swarm % Particle(k) % proc = 0
      Swarm % Particle(k) % buff = 0

      Swarm % Particle(k) % x_o = Swarm % Particle(k) % x_n
      Swarm % Particle(k) % y_o = Swarm % Particle(k) % y_n
      Swarm % Particle(k) % z_o = Swarm % Particle(k) % z_n

      ! Searching for the closest cell and node to place the moved Particle
      call Swarm % Particle(k) % Find_Nearest_Cell(n_parts_in_buffers)
      call Swarm % Particle(k) % Find_Nearest_Node()

      c = Swarm % Particle(k) % cell

      ! Set initial Particle velocities
      rx = Swarm % Particle(k) % x_n - Grid % xc(c)
      ry = Swarm % Particle(k) % y_n - Grid % yc(c)
      rz = Swarm % Particle(k) % z_n - Grid % zc(c)

      ! Compute velocities at the Particle position from velocity gradients
      Swarm % Particle(k) % u    &
         = Flow % u % n(c)       &  ! u velocity at the new time step (% n)
         + Flow % u % x(c) * rx  &  ! u % x is gradient du/dx
         + Flow % u % y(c) * ry  &  ! u % y is gradient du/dy
         + Flow % u % z(c) * rz     ! u % x is gradient du/dz

      Swarm % Particle(k) % v    &
         = Flow % v % n(c)       &  ! v velocity at the new time step (% n)
         + Flow % v % x(c) * rx  &  ! v % x is gradient dv/dx
         + Flow % v % y(c) * ry  &  ! v % y is gradient dv/dy
         + Flow % v % z(c) * rz     ! v % x is gradient dv/dz

      Swarm % Particle(k) % w    &
         = Flow % w % n(c)       &  ! w velocity at the new time step (% n)
         + Flow % w % x(c) * rx  &  ! w % x is gradient dw/dx
         + Flow % w % y(c) * ry  &  ! w % y is gradient dw/dy
         + Flow % w % z(c) * rz     ! w % x is gradient dw/dz

    end do

    Swarm % n_particles = n_new

  end if

  end subroutine
