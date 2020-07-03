!==============================================================================!
  subroutine Swarm_Mod_Allocate(swarm, turb)
!------------------------------------------------------------------------------!
!   Allocates memory to store the charge of each particle                      !
!   It assumes that the number of particles was read from the control file     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
  type(Turb_Type),  target :: turb
!----------------------------------[Locals]------------------------------------!
  integer :: k, nb, nc
!==============================================================================!

  ! Take aliases to object particle flow around
  swarm % pnt_turb => turb
  swarm % pnt_flow => turb % pnt_flow
  swarm % pnt_grid => turb % pnt_grid

  ! Allocate memory for all of them
  if(swarm % n_particles > 0) then
    allocate(swarm % particle(swarm % n_particles))
  end if

  ! Allocate logical array if cell holds particles
  allocate(swarm % cell_has_particles(swarm % pnt_grid % n_cells))
  swarm % cell_has_particles(:) = .false.

  ! Allocate memory for working arrays
  allocate(i_work(swarm % n_particles * N_I_VARS))
  allocate(l_work(swarm % n_particles * N_L_VARS))
  allocate(r_work(swarm % n_particles * N_R_VARS))

  !------------------------------!
  !   Initialize all particles   !
  !------------------------------!
  do k = 1, swarm % n_particles

    ! Take diameter and density from the swarm
    swarm % particle(k) % d       = swarm % diameter
    swarm % particle(k) % density = swarm % density

    ! Set initial velocity to zero
    swarm % particle(k) % u = 0.0
    swarm % particle(k) % v = 0.0
    swarm % particle(k) % w = 0.0

    ! Set initial mean velocity to zero
    swarm % particle(k) % u_mean = 0.0
    swarm % particle(k) % v_mean = 0.0
    swarm % particle(k) % w_mean = 0.0

    ! Set relative velocities to zero (DRW model)
    swarm % particle(k) % rel_u_mod = 0.0
    swarm % particle(k) % rel_v_mod = 0.0
    swarm % particle(k) % rel_w_mod = 0.0

    ! Set DRW velocities to zero (produced by SEIM and seen by particle) 
    swarm % particle(k) % u_drw = 0.0
    swarm % particle(k) % v_drw = 0.0
    swarm % particle(k) % w_drw = 0.0

    ! Set initial coordinates to zero
    swarm % particle(k) % x_n = 0.0
    swarm % particle(k) % y_n = 0.0
    swarm % particle(k) % z_n = 0.0

    swarm % particle(k) % x_o = 0.0
    swarm % particle(k) % y_o = 0.0
    swarm % particle(k) % z_o = 0.0

    ! Set initial cell, node and boundary cell to zero
    swarm % particle(k) % cell     = 0
    swarm % particle(k) % node     = 0
    swarm % particle(k) % bnd_cell = 0
    swarm % particle(k) % bnd_face = 0

    ! Assume particle is in the domain
    ! (A smarter way could be worked out, depending ...
    ! ... on the result of the call to Find_Nearest_Cell)
    swarm % particle(k) % deposited = .false.
    swarm % particle(k) % escaped   = .false.

    ! Is particle in this processor?
    swarm % particle(k) % proc = 0
    swarm % particle(k) % buff = 0

    ! Set initial number of particle states to zero
    swarm % particle(k) % n_states = 0

  end do

  ! Allocate variables for ensemble-averaging
  nb = turb % pnt_grid % n_bnd_cells
  nc = turb % pnt_grid % n_cells
  allocate(swarm % u_mean  (-nb:nc));  swarm % u_mean   = 0.
  allocate(swarm % v_mean  (-nb:nc));  swarm % v_mean   = 0.
  allocate(swarm % w_mean  (-nb:nc));  swarm % w_mean   = 0.
  allocate(swarm % uu      (-nb:nc));  swarm % uu       = 0.
  allocate(swarm % vv      (-nb:nc));  swarm % vv       = 0.
  allocate(swarm % ww      (-nb:nc));  swarm % ww       = 0.
  allocate(swarm % uv      (-nb:nc));  swarm % uv       = 0.
  allocate(swarm % uw      (-nb:nc));  swarm % uw       = 0.
  allocate(swarm % vw      (-nb:nc));  swarm % vw       = 0.
  allocate(swarm % n_states(-nb:nc));  swarm % n_states = 0

  ! Allocate Brownnian diffusion force components
  allocate(swarm % f_fuka_x(-nb:nc));  swarm % f_fuka_x = 0.
  allocate(swarm % f_fuka_y(-nb:nc));  swarm % f_fuka_y = 0.
  allocate(swarm % f_fuka_z(-nb:nc));  swarm % f_fuka_z = 0.

  ! Allocate variables for the modeled flow quantity "v^2"
  allocate(swarm % w_mod_x(-nb:nc));  swarm % w_mod_x   = 0.
  allocate(swarm % w_mod_y(-nb:nc));  swarm % w_mod_y   = 0.
  allocate(swarm % w_mod_z(-nb:nc));  swarm % w_mod_z   = 0.

  ! Gradients of wall-normal modele velocity (retrieved from zeta)
  allocate(swarm % w_mod(-nb:nc));  swarm % w_mod = 0.
  allocate(swarm % w_x  (-nb:nc));  swarm % w_x   = 0.
  allocate(swarm % w_y  (-nb:nc));  swarm % w_y   = 0.
  allocate(swarm % w_z  (-nb:nc));  swarm % w_z   = 0.

  end subroutine
