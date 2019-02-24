!==============================================================================!
  subroutine Swarm_Mod_Create(flow, swarm)
!------------------------------------------------------------------------------!
!   Allocates memory to store the charge of each particle                      !
!   It assumes that the number of particles was read from the control file     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  type(Swarm_Type), target :: swarm
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
    swarm % particle(k) % x = 0.0
    swarm % particle(k) % y = 0.0
    swarm % particle(k) % z = 0.0

    ! Set initial cell, node and boundary cell to zero
    swarm % particle(k) % cell     = 0
    swarm % particle(k) % node     = 0
    swarm % particle(k) % bnd_cell = 0

    ! Assume particle is in the domain
    ! (A smarter way could be worked out, depending ...
    ! ... on the result of the call to Find_Nearest_Cell)
    swarm % particle(k) % deposited = .false.
    swarm % particle(k) % escaped   = .false.

  end do

  end subroutine
