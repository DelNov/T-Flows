!==============================================================================!
  subroutine Initialize_Particle(Particle, Flow, diameter, density)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Particle_Type)     :: Particle
  type(Field_Type), target :: Flow
  real                     :: diameter
  real                     :: density
!==============================================================================!

  ! Store the grid pointer
  Particle % pnt_flow => Flow

  ! Call parent's constructor
  call Particle % Initialize_Point(Flow % pnt_grid)

  ! Take diameter and density from the sender (swarm, most likely)
  Particle % d       = diameter
  Particle % density = density

  ! Set initial velocity to zero
  Particle % u = 0.0
  Particle % v = 0.0
  Particle % w = 0.0

  ! Set relative velocities to zero (DRW model)
  Particle % rel_u_mod = 0.0
  Particle % rel_v_mod = 0.0
  Particle % rel_w_mod = 0.0

  ! Set DRW velocities to zero (produced by SEIM and seen by particle)
  Particle % u_drw = 0.0
  Particle % v_drw = 0.0
  Particle % w_drw = 0.0

  ! Set initial coordinates to zero
  Particle % x_n = 0.0
  Particle % y_n = 0.0
  Particle % z_n = 0.0

  Particle % x_o = 0.0
  Particle % y_o = 0.0
  Particle % z_o = 0.0

  ! Set initial cell, node and boundary cell to zero
  Particle % cell     = 0
  Particle % node     = 0
  Particle % bnd_cell = 0

  ! Assume particle is in the domain
  ! (A smarter way could be worked out, depending ...
  ! ... on the result of the call to Find_Nearest_Cell)
  Particle % deposited = .false.
  Particle % escaped   = .false.

  ! Set some processor number to particle
  Particle % proc = min(1, n_proc)
  Particle % buff = min(1, n_proc)

  end subroutine
