!==============================================================================!
  subroutine Swarm_Mod_Advance_Particles(swarm, turb)
!------------------------------------------------------------------------------!
!   Advances all particles in the swarm.                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
  type(Turb_Type),  target :: turb
  integer                  :: k      ! particle number
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: grid
  type(Field_Type),    pointer :: flow
  type(Particle_Type), pointer :: part
  logical,             pointer :: escaped
  logical,             pointer :: deposited
  integer                      :: ss         ! sub-step counter
!==============================================================================!

  ! Take aliases for the swarm
  grid => swarm % pnt_grid
  flow => swarm % pnt_flow

  ! Take aliases for the particle
  part      => swarm % particle(k)
  escaped   => part  % escaped
  deposited => part  % deposited

  ! Particle time step (division of the global time step)
  swarm % dt = flow % dt / swarm % n_sub_steps

  !----------------------------------!
  !                                  !
  !   Browse through all particles   !
  !                                  !
  !----------------------------------!
  do ss = 1, swarm % n_sub_steps
    do k = 1, swarm % n_particles

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

          ! Make sure particle didn't run out of periodicity
          call Swarm_Mod_Check_Periodicity(swarm, k)
          if(part % node .eq. 0) then
            call Swarm_Mod_Find_Nearest_Cell(swarm, k)
            call Swarm_Mod_Find_Nearest_Node(swarm, k)
          end if

          ! Compute velocity at the particle, and move it
          ! (also calls Bounce_Particle)
          call Swarm_Mod_Move_Particle(swarm, turb, k)

          ! Calling particle forces subroutine to ...
          ! ... compute the forces on each particle and store it
          call Swarm_Mod_Particle_Forces(swarm, k)

        end if  ! in this processor
      end if    ! not deposited or escaped
    end do      ! through particles
  end do        ! through sub-steps

  !---------------------------------------------!
  !                                             !
  !   Exchange particles for parallel version   !
  !                                             !
  !---------------------------------------------!
  call Swarm_Mod_Exchange_Particles(swarm)

  !-----------------------------------!
  !   Print some data on the screen   !
  !-----------------------------------!
  if(this_proc < 2) then
    do k = 1, swarm % n_particles

      ! Printing particle position
      write(*,'(a,i3,a,i7,a,i2,a,3es15.6,a,es12.4)')                 &
              '#  particle: ',  k,                                   &
              ',  cell: ',      part % cell,                         &
              ',  processor: ', part % proc,                         &
              ',  x,y,z: ',     part % x_n, part % y_n, part % z_n,  &
              ',  cfl: ',       part % cfl
    end do
  end if

  end subroutine
