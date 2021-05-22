!==============================================================================!
  subroutine Swarm_Mod_Allocate(swarm, Flow, turb)
!------------------------------------------------------------------------------!
!   Allocates memory to store the charge of each Particle                      !
!   It assumes that the number of particles was read from the control file     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: turb
!----------------------------------[Locals]------------------------------------!
  integer :: k, nb, nc
!==============================================================================!

  ! Take aliases to object Particle Flow around
  swarm % pnt_flow => Flow
  swarm % pnt_grid => Flow % pnt_grid
  swarm % pnt_turb => turb

  ! Allocate memory for all of them
  if(swarm % n_particles > 0) then
    allocate(swarm % Particle(swarm % n_particles))
  end if

  ! Allocate logical array if cell holds particles
  allocate(swarm % cell_has_particles(swarm % pnt_grid % n_cells))
  swarm % cell_has_particles(:) = .false.

  ! Allocate memory for working arrays
  allocate(swarm % i_work(swarm % n_particles * swarm % N_I_VARS))
  allocate(swarm % l_work(swarm % n_particles * swarm % N_L_VARS))
  allocate(swarm % r_work(swarm % n_particles * swarm % N_R_VARS))

  !------------------------------!
  !   Initialize all particles   !
  !------------------------------!
  do k = 1, swarm % n_particles
    call swarm % Particle(k) % Initialize_Particle(Flow, swarm % diameter,  &
                                                         swarm % density)
  end do

  ! Aliases for cell-based variables
  nb = turb % pnt_grid % n_bnd_cells
  nc = turb % pnt_grid % n_cells

  ! Reflected and deposited particles on the walls and the escaped particles
  if(nb > 0) then
    allocate(swarm % n_reflected(-nb:nc));  swarm % n_reflected(:) = 0
    allocate(swarm % n_deposited(-nb:nc));  swarm % n_deposited(:) = 0
    allocate(swarm % n_escaped  (-nb:nc));  swarm % n_escaped(:)   = 0
  else  ! take care of subdomains which have no boundary cells
    allocate(swarm % n_reflected(0:0));  swarm % n_reflected(:) = 0
    allocate(swarm % n_deposited(0:0));  swarm % n_deposited(:) = 0
    allocate(swarm % n_escaped  (0:0));  swarm % n_escaped(:)   = 0
  end if

  ! Allocate variables for ensemble-averaging
  if(swarm % statistics) then
    allocate(swarm % u_mean  (-nb:nc));  swarm % u_mean(:)   = 0.
    allocate(swarm % v_mean  (-nb:nc));  swarm % v_mean(:)   = 0.
    allocate(swarm % w_mean  (-nb:nc));  swarm % w_mean(:)   = 0.
    allocate(swarm % uu      (-nb:nc));  swarm % uu(:)       = 0.
    allocate(swarm % vv      (-nb:nc));  swarm % vv(:)       = 0.
    allocate(swarm % ww      (-nb:nc));  swarm % ww(:)       = 0.
    allocate(swarm % uv      (-nb:nc));  swarm % uv(:)       = 0.
    allocate(swarm % uw      (-nb:nc));  swarm % uw(:)       = 0.
    allocate(swarm % vw      (-nb:nc));  swarm % vw(:)       = 0.
    allocate(swarm % n_states(-nb:nc));  swarm % n_states(:) = 0
  end if

  ! Allocate Brownnian diffusion force components
  allocate(swarm % f_fuka_x(-nb:nc));  swarm % f_fuka_x(:) = 0.
  allocate(swarm % f_fuka_y(-nb:nc));  swarm % f_fuka_y(:) = 0.
  allocate(swarm % f_fuka_z(-nb:nc));  swarm % f_fuka_z(:) = 0.

  ! Allocate variables for the modeled Flow quantity "v^2"
  allocate(swarm % v2_mod  (-nb:nc));  swarm % v2_mod(:)   = 0.
  allocate(swarm % v2_mod_x(-nb:nc));  swarm % v2_mod_x(:) = 0.
  allocate(swarm % v2_mod_y(-nb:nc));  swarm % v2_mod_y(:) = 0.
  allocate(swarm % v2_mod_z(-nb:nc));  swarm % v2_mod_z(:) = 0.

  end subroutine
