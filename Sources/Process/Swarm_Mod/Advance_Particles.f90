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
  real                         :: pr(4)      ! stuff to print
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

      ! Calling the nearest cell subroutine to find the ...
      ! ... nearest cell for each particle and stores it
      call Swarm_Mod_Find_Nearest_Cell(swarm, k)

      ! If particle is in this processor, carry on with it
      if(part % here) then

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

      ! Print particle position and cfl number
      pr(1:4) = 0.0
      if(part % here) then
        pr(1) = part % x_n
        pr(2) = part % y_n
        pr(3) = part % z_n
        pr(4) = part % cfl
      end if
      call Comm_Mod_Global_Sum_Real_Array(pr, 4)
      if(this_proc < 2) then
        ! Printing particle position
        print *,k, 'position','(', pr(1), ',',  &
                                   pr(2), ',',  &
                                   pr(3), ')',  &
                   ', | cfl =',    pr(4)
      end if

    end if  ! not deposited and not escaped

  end do

  end subroutine
