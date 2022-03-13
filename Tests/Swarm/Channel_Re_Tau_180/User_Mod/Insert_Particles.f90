!==============================================================================!
  subroutine User_Mod_Insert_Particles(Flow, Turb, Vof, Swarm, n, time)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of time step.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
  integer, intent(in)      :: n     ! time step
  real,    intent(in)      :: time  ! physical time
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: c, i, j, k, n_parts_in_buffers
  real                     :: x, y, z, dy, dz, my, mz
  real                     :: rx, ry, rz
  real                     :: c1, c2, c3   ! random variables
!------------------------------[Local parameters]------------------------------!
  real, parameter :: L1 = 6.28  ! streamwise
  real, parameter :: L2 = 3.14  ! spanwise
  real, parameter :: L3 = 2.0   ! wall-normal
!==============================================================================!

  ! Take alias(es)
  Grid => Flow % pnt_grid

  ! Random variables
  c1 = 1.0
  c2 = 1.0
  c3 = 1.0

  !-------------------!
  !   1st time step   !
  !-------------------!
  if(n .eq. 1) then     ! should be after the Flow is developed

    ! Browsing through all introduced particles
    do k = 1, Swarm % n_particles

        ! Initalizing particle position (already initialized in
        ! Swarm_Mod_Allocate)
        Swarm % Particle(k) % x_n = 0.0
        Swarm % Particle(k) % y_n = 0.0
        Swarm % Particle(k) % z_n = 0.0

        ! Generating random locations for particle
        call random_number(c1)
        call random_number(c2)
        call random_number(c3)

        ! Initalizing Particle position
        Swarm % Particle(k) % x_n = (L1 * c1) + Swarm % Particle(k) % x_n
        Swarm % Particle(k) % y_n = (L2 * c2) + Swarm % Particle(k) % y_n
        Swarm % Particle(k) % z_n = (L3 * c3) + Swarm % Particle(k) % z_n

        ! you essentially moved them a lot (from 0, 0, 0)
        Swarm % Particle(k) % cell = 0
        Swarm % Particle(k) % node = 0
        Swarm % Particle(k) % proc = 0
        Swarm % Particle(k) % buff = 0

        Swarm % Particle(k) % x_o = Swarm % Particle(k) % x_n
        Swarm % Particle(k) % y_o = Swarm % Particle(k) % y_n
        Swarm % Particle(k) % z_o = Swarm % Particle(k) % z_n

        ! Searching for the closest cell and node to place the moved particle
        call Swarm % Particle(k) % Find_Nearest_Cell(n_parts_in_buffers)
        call Swarm % Particle(k) % Find_Nearest_Node()

        c = Swarm % Particle(k) % cell

        ! Set initial Particle velocities
        rx = Swarm % Particle(k) % x_n - Grid % xc(c)
        ry = Swarm % Particle(k) % y_n - Grid % yc(c)
        rz = Swarm % Particle(k) % z_n - Grid % zc(c)

        ! Compute velocities at the particle position from velocity gradients
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
  end if

  end subroutine
