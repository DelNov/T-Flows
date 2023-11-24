!==============================================================================!
  subroutine Insert_At(Particle, x, y, z, n_parts_in_buffers, Flow, Vof)
!------------------------------------------------------------------------------!
!> The Insert_At subroutine repositions a particle to a specified location
!> within the computational grid and updates its associated properties. It is
!> crucial for accurately placing particles and initializing their state in
!> simulations involving Lagrangian particle tracking within the Swarm_Mod.
!------------------------------------------------------------------------------!
! Functionality                                                                !
!                                                                              !
! * Particle repositioning: Updates the particle's position to the specified   !
!   coordinates, essential for dynamic particle simulations.                   !
! * Closest cell and node determination: Finds the nearest cell and node to    !
!   the new position of the particle, ensuring accurate spatial referencing.   !
! * Velocity computation: Calculates particle's velocity based on the flow     !
!   field at its new position, essential for fluid-particle interaction.       !
! * Smoothed VOF value calculation: Determines the smoothed Volume of Fluid    !
!   (VOF) value at the particle's position, aiding in fluid interaction        !
!   modeling for three-phase flow simulations.                                 !
! * Previous position update: Updates the particle's old position to match     !
!   the new one, maintaining consistency in tracking particle's trajectory.    !
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Particle_Type)       :: Particle            !! particle being inserted
  real,    intent(in)        :: x, y, z             !! particle coordinate
  integer                    :: n_parts_in_buffers  !! counter for the number of
                                                    !! particles in buffer zones
  type(Field_Type), optional :: Flow                !! flow in which particle is
                                                    !! being inserted
  type(Vof_Type),   optional :: Vof                 !! VOF for three-phase flows
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: Grid
  real                     :: rx, ry, rz
  integer                  :: c
!==============================================================================!

  ! Store the grid pointer
  Grid => Particle % pnt_grid

  Particle % x_n = x
  Particle % y_n = y
  Particle % z_n = z

  ! You essentially moved them a lot (from 0, 0, 0)
  Particle % cell = 0
  Particle % node = 0
  Particle % proc = 0
  Particle % buff = 0

  ! Searching for the closest cell and node to place the moved particle
  call Particle % Find_Nearest_Cell(n_parts_in_buffers)
  call Particle % Find_Nearest_Node()

  ! Take an alias
  c = Particle % cell

  ! Set initial particle velocities
  rx = Particle % x_n - Grid % xc(c)
  ry = Particle % y_n - Grid % yc(c)
  rz = Particle % z_n - Grid % zc(c)

  ! Compute velocities at the particle position from velocity gradients
  if(present(Flow)) then
    Particle % u               &
       = Flow % u % n(c)       &  ! u velocity at the new time step (% n)
       + Flow % u % x(c) * rx  &  ! u % x is gradient du/dx
       + Flow % u % y(c) * ry  &  ! u % y is gradient du/dy
       + Flow % u % z(c) * rz     ! u % x is gradient du/dz

    Particle % v               &
       = Flow % v % n(c)       &  ! v velocity at the new time step (% n)
       + Flow % v % x(c) * rx  &  ! v % x is gradient dv/dx
       + Flow % v % y(c) * ry  &  ! v % y is gradient dv/dy
       + Flow % v % z(c) * rz     ! v % x is gradient dv/dz

    Particle % w               &
       = Flow % w % n(c)       &  ! w velocity at the new time step (% n)
       + Flow % w % x(c) * rx  &  ! w % x is gradient dw/dx
       + Flow % w % y(c) * ry  &  ! w % y is gradient dw/dy
       + Flow % w % z(c) * rz     ! w % x is gradient dw/dz
  end if

  ! Value of smoothed vof at the old position
  if(present(Vof)) then
    Particle % smooth_n = Vof % smooth % n(c)       &
                        + Vof % smooth % x(c) * rx  &
                        + Vof % smooth % y(c) * ry  &
                        + Vof % smooth % z(c) * rz
    Particle % smooth_o = Particle % smooth_n
  end if

  ! Or move back against velocity?
  Particle % x_o = Particle % x_n
  Particle % y_o = Particle % y_n
  Particle % z_o = Particle % z_n

  end subroutine
