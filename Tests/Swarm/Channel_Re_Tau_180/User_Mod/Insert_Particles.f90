!==============================================================================!
  subroutine User_Mod_Insert_Particles(flow, turb, mult, swarm, n, time)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of time step.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Swarm_Type),      target :: swarm
  integer                       :: n     ! time step
  real                          :: time  ! physical time
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: c, i, j, k, n_parts_in_buffers
  real                     :: x, y, z, dy, dz, my, mz
  real                     :: rx, ry, rz
  real                     :: c1, c2, c3   ! random variables
!------------------------------[Local parameters]------------------------------!
  real, parameter :: L1 = 6.28  ! streamwise
  real, parameter :: L2 = 3.14  ! spanwise
  real, parameter :: L3 = 2.0   ! wall-normal
!==============================================================================!

  ! Random variables
  c1 = 1.0
  c2 = 1.0
  c3 = 1.0

  !-------------------!
  !   1st time step   !
  !-------------------!
  if(n .eq. 1) then     ! should be after the flow is developed

    ! Browsing through all introduced particles
    do k = 1, swarm % n_particles

        ! Initalizing particle position (already initialized in
        ! Swarm_Mod_Allocate)
        swarm % particle(k) % x_n = 0.0
        swarm % particle(k) % y_n = 0.0
        swarm % particle(k) % z_n = 0.0

        ! Generating random locations for particle
        call random_number(c1)
        call random_number(c2)
        call random_number(c3)

        ! Initalizing particle position
        swarm % particle(k) % x_n = (L1 * c1) + swarm % particle(k) % x_n
        swarm % particle(k) % y_n = (L2 * c2) + swarm % particle(k) % y_n
        swarm % particle(k) % z_n = (L3 * c3) + swarm % particle(k) % z_n

        ! you essentially moved them a lot (from 0, 0, 0)
        swarm % particle(k) % cell = 0
        swarm % particle(k) % node = 0
        swarm % particle(k) % proc = 0
        swarm % particle(k) % buff = 0

        swarm % particle(k) % x_o = swarm % particle(k) % x_n
        swarm % particle(k) % y_o = swarm % particle(k) % y_n
        swarm % particle(k) % z_o = swarm % particle(k) % z_n

        ! Searching for the closest cell and node to place the moved particle
        call Swarm_Mod_Find_Nearest_Cell(swarm, k, n_parts_in_buffers)
        call Swarm_Mod_Find_Nearest_Node(swarm, k)

        c = swarm % particle(k) % cell

        ! Set initial particle velocities
        rx = swarm % particle(k) % x_n - grid % xc(c)
        ry = swarm % particle(k) % y_n - grid % yc(c)
        rz = swarm % particle(k) % z_n - grid % zc(c)

        ! Compute velocities at the particle position from velocity gradients
        swarm % particle(k) % u    &
           = flow % u % n(c)       &  ! u velocity at the new time step (% n)
           + flow % u % x(c) * rx  &  ! u % x is gradient du/dx
           + flow % u % y(c) * ry  &  ! u % y is gradient du/dy
           + flow % u % z(c) * rz     ! u % x is gradient du/dz

        swarm % particle(k) % v    &
           = flow % v % n(c)       &  ! v velocity at the new time step (% n)
           + flow % v % x(c) * rx  &  ! v % x is gradient dv/dx
           + flow % v % y(c) * ry  &  ! v % y is gradient dv/dy
           + flow % v % z(c) * rz     ! v % x is gradient dv/dz

        swarm % particle(k) % w    &
           = flow % w % n(c)       &  ! w velocity at the new time step (% n)
           + flow % w % x(c) * rx  &  ! w % x is gradient dw/dx
           + flow % w % y(c) * ry  &  ! w % y is gradient dw/dy
           + flow % w % z(c) * rz     ! w % x is gradient dw/dz

    end do
  end if

  end subroutine
