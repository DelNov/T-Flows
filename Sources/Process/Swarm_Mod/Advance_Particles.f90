!==============================================================================!
  subroutine Swarm_Mod_Advance_Particles(swarm)
!------------------------------------------------------------------------------!
!   Advances all particles in the swarm.                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
  integer                  :: k      ! particle number
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: grid
  type(Particle_Type), pointer :: part
  logical,             pointer :: escaped
  logical,             pointer :: deposited
!==============================================================================!

  ! Take aliases for the swarm
  grid => swarm % pnt_grid

  !----------------------------------!
  !   Browse through all particles   !
  !----------------------------------!
  do k = 1, swarm % n_particles

    ! Take aliases for the particle
    part      => swarm % particle(k)
    escaped   => part  % escaped
    deposited => part  % deposited

    ! If particle is neither deposited nor escped
    if(.not. deposited .and. .not. escaped) then

      ! Calling the nearest cell subroutine to find the ...
      ! ... nearest cell for each particle and stores it
      call Swarm_Mod_Find_Nearest_Cell(swarm, k)

      ! Calling the nearest node subroutine to find the ...
      ! ... nearest node for each particle and stores it
      call Swarm_Mod_Find_Nearest_Node(swarm, k)

      ! Calling interpolate velocity subroutine ...
      ! ... to update velocity of each particle
      call Swarm_Mod_Move_Particle(swarm, k)

      ! Calling particle forces subroutine to ...
      ! ... compute the forces on each particle and store it
      call Swarm_Mod_Particle_Forces(swarm, k)
    end if

  end do

  end subroutine
