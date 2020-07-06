!==============================================================================!
  subroutine Swarm_Mod_Advance_Particles(swarm, n, n_stat_p, first_dt_p)
!------------------------------------------------------------------------------!
!   Advances all particles in the swarm.                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
  integer, intent(in)      :: n           ! current time step
  integer, intent(in)      :: n_stat_p    ! starting time for swarm statistics
  integer, intent(in)      :: first_dt_p  ! starting time for swarm simulation
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: grid
  type(Field_Type),    pointer :: flow
  type(Turb_Type),     pointer :: turb
  type(Particle_Type), pointer :: part
  logical,             pointer :: escaped
  logical,             pointer :: deposited
  integer                      :: k                   ! particle number
  integer                      :: ss                  ! sub-step counter
  integer                      :: n_parts_in_buffers
  real                         :: avg_part_cfl, avg_part_re, avg_part_st
  real                         :: max_part_cfl, max_part_re, max_part_st
!==============================================================================!

  ! Take aliases for the swarm
  grid => swarm % pnt_grid
  flow => swarm % pnt_flow
  turb => swarm % pnt_turb

  ! Particle time step (division of the global time step)
  swarm % dt = flow % dt / swarm % n_sub_steps

  !------------------------!
  !   Fukagata SGS model   !
  !------------------------!
  if(turb % model .eq. HYBRID_LES_PRANDTL) then
    if(swarm % subgrid_scale_model .eq. BROWNIAN_FUKAGATA) then
      call Swarm_Mod_Sgs_Fukagata(swarm)
    end if
  end if

  if(turb % model .eq. HYBRID_LES_RANS) then

    ! Correcting for particle time step size (if ER-HRL model is used)
    call Swarm_Mod_Particle_Time_Scale(swarm)

    ! Store gradients for modeled flow quantities for swarm
    call Swarm_Mod_Grad_Modeled_Flow(swarm, k)

  end if

  ! Gaussian random no.s interval (for SEIM model)
  swarm % time_eim = n - first_dt_p

  !---------------------------------------------!
  !       Store old particle coordinates        !
  !   (important for bouncing from the walls)   !
  !---------------------------------------------!
  do k = 1, swarm % n_particles
    part => swarm % particle(k)
    part % x_o = part % x_n
    part % y_o = part % y_n
    part % z_o = part % z_n
  end do

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

      !--------------------------------------------------!
      !   If particle is neither deposited nor escaped   !
      !--------------------------------------------------!
      if(.not. deposited .and. .not. escaped) then

        ! If particle is in this processor, carry on with it
        if(part % proc .eq. this_proc) then

          ! Compute velocity at the particle, and move it
          ! (also calls Bounce_Particle)
          call Swarm_Mod_Move_Particle(swarm, k)

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
          call Swarm_Mod_Calculate_Mean(swarm, k, n, n_stat_p, ss)

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
  call Swarm_Mod_Print_Statistics(swarm)

  end subroutine
