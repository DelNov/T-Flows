!==============================================================================!
  subroutine Insert_At(Particle, x, y, z, n_parts_in_buffers)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Particle_Type) :: Particle
  real,    intent(in)  :: x, y, z
  integer              :: n_parts_in_buffers
!----------------------------------[Locals]------------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  real                      :: rx, ry, rz
  integer                   :: c
!==============================================================================!

  ! Store the grid pointer
  Flow => Particle % pnt_flow
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

  ! Or move back against velocity?
  Particle % x_o = Particle % x_n
  Particle % y_o = Particle % y_n
  Particle % z_o = Particle % z_n

  end subroutine
