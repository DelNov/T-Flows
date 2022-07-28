!==============================================================================!
  subroutine Allocate_Swarm(Swarm, Flow, Turb, Vof)
!------------------------------------------------------------------------------!
!   Allocates memory to store the charge of each particle                      !
!   It assumes that the number of particles was read from the control file     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Swarm_Type), target :: Swarm
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: Turb
  type(Vof_Type),    target :: Vof
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
    if(this_proc < 2) then
      print *, '# ERROR!  You are attempting a simulation with'
      print *, '# particles but max number of particles is zero.'
      print *, '# Did you set the parameter MAX_PARTICLES in control file?'
      print *, '# This error is critical.  Exiting!'
      call Comm_Mod_End
      stop
    end if
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
