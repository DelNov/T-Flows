!==============================================================================!
  subroutine Swarm_Mod_Advance_Particles(swarm, turb, n, n_stat_p)
!------------------------------------------------------------------------------!
!   Advances all particles in the swarm.                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
  type(Turb_Type),  target :: turb
  integer                  :: n          ! current time step
  integer                  :: n_stat_p   ! starting time for swarm statistics
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: grid
  type(Field_Type),    pointer :: flow
  type(Particle_Type), pointer :: part
  logical,             pointer :: escaped
  logical,             pointer :: deposited
  integer                      :: k          ! particle number
  integer                      :: ss         ! sub-step counter
  integer                      :: n_parts_in_buffers
!==============================================================================!

  ! Take aliases for the swarm
  grid => swarm % pnt_grid
  flow => swarm % pnt_flow

  ! Particle time step (division of the global time step)
  swarm % dt = flow % dt / swarm % n_sub_steps

  !----------------------------------!
  !                                  !
  !   Browse through all particles   !
  !                                  !
  !----------------------------------!
  do ss = 1, swarm % n_sub_steps

    n_parts_in_buffers = 0

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

          ! Compute velocity at the particle, and move it
          ! (also calls Bounce_Particle)
          call Swarm_Mod_Move_Particle(swarm, turb, k)

          ! Calling particle forces subroutine to ...
          ! ... compute the forces on each particle and store it
          call Swarm_Mod_Particle_Forces(swarm, k)

          ! Calling the nearest cell subroutine to find the ...
          ! ... nearest cell for each particle and stores it
          call Swarm_Mod_Find_Nearest_Cell(swarm, k, n_parts_in_buffers)

          ! Calling the nearest node subroutine to find the ...
          ! ... nearest node for each particle and stores it
          call Swarm_Mod_Find_Nearest_Node(swarm, k)

          ! First check if it didn't escape through periodicity
          call Swarm_Mod_Check_Periodicity(swarm, k, n_parts_in_buffers)

          ! Gathering swarm statistics  
          call Swarm_Mod_Calculate_Mean(swarm, k, n, n_stat_p)

        end if  ! in this processor
      end if    ! deposited or escaped
    end do      ! through particles

    ! Exchange particles for parallel version; if needed
    call Comm_Mod_Global_Sum_Int(n_parts_in_buffers)
    if(n_parts_in_buffers > 0) then
      call Swarm_Mod_Exchange_Particles(swarm)
    end if

  end do        ! through sub-steps


  !-----------------------------------!
  !   Print some data on the screen   !
  !-----------------------------------!
  do k = 1, swarm % n_particles

    ! Refresh the alias
    part => swarm % particle(k)
    
    if(this_proc .eq. part % proc) then
      ! Printing particle position
      write(*,'(a,i7,a,i7,a,i2,a,3es15.6,a,es12.4)')                 &
              '# particle:',  k,                                     &
              ',  cell: ',      part % cell,                         &
              ',  processor: ', part % proc,                         &
              ',  x,y,z: ',     part % x_n, part % y_n, part % z_n,  &
              ',  cfl: ',       part % cfl
    end if
  end do

  end subroutine
