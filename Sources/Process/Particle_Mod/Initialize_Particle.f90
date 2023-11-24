!==============================================================================!
  subroutine Initialize_Particle(Particle, Flow, diameter, density)
!------------------------------------------------------------------------------!
!> The Initialize_Particle subroutine is tasked with setting up the initial
!> state of a particle in the computational domain. It initializes the
!> particle's properties and links it to the flow field, ensuring that the
!> particle is ready for simulation and interaction within the fluid
!> environment.
!------------------------------------------------------------------------------!
! Functionality                                                                !
!                                                                              !
! * Linking to flow field: Associates the particle with a specific flow        !
!   field, enabling interaction and tracking within that environment.          !
! * Particle properties: Sets the particle's diameter and density based on     !
!   provided values, crucial for physical simulations and force calculations.  !
! * Velocity initialization: Initializes particle's velocity components to     !
!   zero, preparing it for dynamic updates during the simulation process.      !
! * Relative velocity setup: Sets up initial relative velocities for DRW       !
!   model interactions, important for stochastic eddy interaction modeling.    !
! * Position and state initialization: Establishes the initial position        !
!   (coordinates) of the particle and resets its state flags (deposited,       !
!   escaped, trapped), ensuring accurate tracking and interaction handling.    !
! * Processor assignment: Assigns a processor number to the particle,          !
!   critical in parallel processing contexts for managing particle location.   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Particle_Type)     :: Particle  !! particle object to be initialized
  type(Field_Type), target :: Flow      !! flow in which particle is immersed
  real,         intent(in) :: diameter  !! particle's diameter
  real,         intent(in) :: density   !! particle's density
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
  Particle % trapped   = .false.

  ! Set some processor number to particle
  Particle % proc = min(1, N_Procs())
  Particle % buff = min(1, N_Procs())

  end subroutine
