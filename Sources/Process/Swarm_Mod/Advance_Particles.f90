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
  type(Grid_Type),     pointer     :: grid
  type(Particle_Type), pointer     :: part
  logical,             pointer     :: escaped
  logical,             pointer     :: deposited
!==============================================================================!

  ! Take aliases for the swarm
  grid => swarm % pnt_grid

  !----------------------------------!
  !                                  !
  !   Browse through all particles   !
  !                                  !
  !----------------------------------!
  do k = 1, swarm % n_particles

    ! Take aliases for the particle
    part      => swarm % particle(k)
    escaped   => part  % escaped
    deposited => part  % deposited

    !-------------------------------------------------!
    !   If particle is neither deposited nor escped   !
    !-------------------------------------------------!
    if(.not. deposited .and. .not. escaped) then

      ! If particle is in this processor, carry on with it
      if(part % proc .eq. this_proc) then

        ! Calling the nearest cell subroutine to find the ...
        ! ... nearest cell for each particle and stores it
        call Swarm_Mod_Find_Nearest_Cell(swarm, k)

        ! Calling the nearest node subroutine to find the ...
        ! ... nearest node for each particle and stores it
        call Swarm_Mod_Find_Nearest_Node(swarm, k)

        ! Compute velocity at the particle, and move it
        ! (also calls Bounce_Particle)
        call Swarm_Mod_Move_Particle(swarm, k)

        ! Calling particle forces subroutine to ...
        ! ... compute the forces on each particle and store it
        call Swarm_Mod_Particle_Forces(swarm, k)

      end if  ! in this processor
    end if    ! not deposited or escaped
  end do      ! through particles

  call Swarm_Mod_Exchange_Particles(swarm)

  !-----------------------------------!
  !   Print some data on the screen   !
  !-----------------------------------!
  if(this_proc < 2) then
    do k = 1, swarm % n_particles

      ! Take alias
      part => swarm % particle(k)

      ! Printing particle position
      write(*,'(a,i3,a,i2,a,3e15.6,a,e12.4)')                        &
              '#  particle: ',  k,                                   &
              ',  processor: ', part % proc,                         &
              ',  x,y,z: ',     part % x_n, part % y_n, part % z_n,  &
              ',  cfl: ',       part % cfl
    end do
  end if

  end subroutine
