!==============================================================================!
  subroutine Create_Swarm(Swarm, Flow, Turb, Vof)
!------------------------------------------------------------------------------!
!>  This subroutine initializes and allocates memory for a swarm of particles,
!>  setting up the necessary components for Lagrangian particle tracking. It
!>  configures the swarm based on predefined parameters and links it to the
!>  corresponding flow, turbulence, and VOF fields.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initialization: Sets up a swarm with a specified number of particles     !
!     and links it to the relevant flow, turbulence, and VOF fields.           !
!   * Memory allocation: Allocates memory for storing particle properties and  !
!     various arrays for tracking particle interactions with the flow domain   !
!     and for parallel processing.                                             !
!   * Particle initialization: Initializes each particle in the swarm with     !
!     predefined properties like diameter and density.                         !
!   * Ensemble-Averaging Preparation: Allocates memory for variables needed    !
!     for ensemble-averaging in statistical studies.                           !
!   * Brownian Motion and Modeled Flow: Prepares variables for Brownian motion !
!     forces and modeled flow quantity calculations.                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Swarm_Type), target :: Swarm  !! the swarm od particles
  type(Field_Type),  target :: Flow   !! flow field
  type(Turb_Type),   target :: Turb   !! turbulence field
  type(Vof_Type),    target :: Vof    !! volume of fluid
!----------------------------------[Locals]------------------------------------!
  integer :: k, nb, nc
!==============================================================================!

  if( .not. Flow % with_particles) return

  ! Take aliases to object particle Flow around
  Swarm % pnt_flow => Flow
  Swarm % pnt_grid => Flow % pnt_grid
  Swarm % pnt_turb => Turb
  Swarm % pnt_vof  => Vof

  ! Allocate memory for all of them
  if(Swarm % max_particles > 0) then
    allocate(Swarm % Particle(Swarm % max_particles))
  else
    call Message % Error(72,                                               &
             'You are attempting a simulation with particles, but max '//  &
             'number of particles is zero.  Did you set the parameter '//  &
             'MAX_PARTICLES in control file?  This error is critical. '//  &
             ' Exiting!', file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  ! Allocate logical array if cell holds particles
  allocate(Swarm % cell_has_particles(Swarm % pnt_grid % n_cells))
  Swarm % cell_has_particles(:) = .false.

  ! Allocate memory for working arrays
  allocate(Swarm % i_work(Swarm % max_particles * Swarm % N_I_VARS))
  allocate(Swarm % l_work(Swarm % max_particles * Swarm % N_L_VARS))
  allocate(Swarm % r_work(Swarm % max_particles * Swarm % N_R_VARS))

  !------------------------------!
  !   Initialize all particles   !
  !------------------------------!
  do k = 1, Swarm % max_particles
    call Swarm % Particle(k) % Initialize_Particle(Flow, Swarm % diameter,  &
                                                         Swarm % density)
  end do

  ! Aliases for cell-based variables
  nb = Turb % pnt_grid % n_bnd_cells
  nc = Turb % pnt_grid % n_cells

  ! Reflected and deposited particles on the walls and the escaped particles
  allocate(Swarm % n_reflected(-nb:nc));  Swarm % n_reflected(:) = 0
  allocate(Swarm % n_deposited(-nb:nc));  Swarm % n_deposited(:) = 0
  allocate(Swarm % n_escaped  (-nb:nc));  Swarm % n_escaped(:)   = 0
  Swarm % n_trapped = 0

  ! Allocate variables for ensemble-averaging
  if(Swarm % statistics) then
    allocate(Swarm % u_mean  (-nb:nc));  Swarm % u_mean(:)   = 0.
    allocate(Swarm % v_mean  (-nb:nc));  Swarm % v_mean(:)   = 0.
    allocate(Swarm % w_mean  (-nb:nc));  Swarm % w_mean(:)   = 0.
    allocate(Swarm % uu      (-nb:nc));  Swarm % uu(:)       = 0.
    allocate(Swarm % vv      (-nb:nc));  Swarm % vv(:)       = 0.
    allocate(Swarm % ww      (-nb:nc));  Swarm % ww(:)       = 0.
    allocate(Swarm % uv      (-nb:nc));  Swarm % uv(:)       = 0.
    allocate(Swarm % uw      (-nb:nc));  Swarm % uw(:)       = 0.
    allocate(Swarm % vw      (-nb:nc));  Swarm % vw(:)       = 0.
    allocate(Swarm % n_states(-nb:nc));  Swarm % n_states(:) = 0
  end if

  ! Allocate Brownnian diffusion force components
  allocate(Swarm % f_fuka_x(-nb:nc));  Swarm % f_fuka_x(:) = 0.
  allocate(Swarm % f_fuka_y(-nb:nc));  Swarm % f_fuka_y(:) = 0.
  allocate(Swarm % f_fuka_z(-nb:nc));  Swarm % f_fuka_z(:) = 0.

  ! Allocate variables for the modeled Flow quantity "v^2"
  allocate(Swarm % v2_mod  (-nb:nc));  Swarm % v2_mod(:)   = 0.
  allocate(Swarm % v2_mod_x(-nb:nc));  Swarm % v2_mod_x(:) = 0.
  allocate(Swarm % v2_mod_y(-nb:nc));  Swarm % v2_mod_y(:) = 0.
  allocate(Swarm % v2_mod_z(-nb:nc));  Swarm % v2_mod_z(:) = 0.

  end subroutine
