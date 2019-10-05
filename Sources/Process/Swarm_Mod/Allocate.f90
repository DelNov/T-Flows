!==============================================================================!
  subroutine Swarm_Mod_Allocate(swarm, flow)
!------------------------------------------------------------------------------!
!   Allocates memory to store the charge of each particle                      !
!   It assumes that the number of particles was read from the control file     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
  type(Field_Type), target :: flow
!----------------------------------[Locals]------------------------------------!
  integer :: k
!==============================================================================!

  ! Take aliases to object particle flow around
  swarm % pnt_flow => flow
  swarm % pnt_grid => flow % pnt_grid

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

  end do

  end subroutine
