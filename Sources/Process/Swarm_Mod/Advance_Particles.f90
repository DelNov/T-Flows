!==============================================================================!
  subroutine Advance_Particles(Swarm, n_stat_p, first_dt_p)
!------------------------------------------------------------------------------!
!> This subroutine advances all particles in a swarm through a computational
!> domain. It handles particle motion, interactions with the flow field,
!> turbulence effects, and interactions with boundaries and surfaces. The
!> routine is a central part of the Lagrangian particle tracking approach in
!> T-Flows.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Particle motion: Manages the movement of particles within the grid based !
!     on the current flow field and particle properties.                       !
!   * Subgrid scale modeling: Applies turbulence models for particle tracking, !
!     including Brownian motion (Fukagata) and Discrete Random Walk (DRW).     !
!   * Boundary interaction: Processes particle collision and interaction with  !
!     computational domain boundaries, including bouncing and trapping.        !
!   * Surface interaction: Manages the behavior of particles upon contact with !
!     surfaces, such as trapping or deposition.                                !
!   * Parallel processing: Supports particle exchange in parallel computation  !
!     environments, ensuring accurate tracking across processors.              !
!   * Statistical analysis: Calculates and updates particle statistics during  !
!     the simulation.                                                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Swarm_Type), target :: Swarm       !! the swarm of particles
  integer, intent(in)       :: n_stat_p    !! starting time for swarm statistics
  integer, intent(in)       :: first_dt_p  !! starting time for swarm simulation
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: Grid
  type(Field_Type),    pointer :: Flow
  type(Turb_Type),     pointer :: Turb
  type(Particle_Type), pointer :: Part
  logical,             pointer :: deposited
  logical,             pointer :: escaped
  logical,             pointer :: trapped
  integer                      :: k                   ! particle number
  integer                      :: ss                  ! sub-step counter
  integer                      :: n_parts_in_buffers
!==============================================================================!

  ! Take aliases for the Swarm
  Grid => Swarm % pnt_grid
  Flow => Swarm % pnt_flow
  Turb => Swarm % pnt_turb

  ! Particle time step (division of the global time step)
  Swarm % dt = Flow % dt / real(Swarm % n_sub_steps)

  !------------------------!
  !   Fukagata SGS model   !
  !------------------------!
  if(Turb % model .eq. HYBRID_LES_PRANDTL) then
    if(Swarm % subgrid_scale_model .eq. BROWNIAN_FUKAGATA) then
      call Swarm % Sgs_Fukagata()
    end if
  end if

  if(Turb % model .eq. HYBRID_LES_RANS) then

    ! Correcting for particle time step size (if ER-HRL model is used)
    call Swarm % Particle_Time_Scale()

    ! Store gradients for modeled Flow quantities for Swarm
    call Swarm % Grad_Modeled_Flow()

  end if

  ! Gaussian random no.s interval (for SEIM model)
  Swarm % time_eim = Time % Curr_Dt() - first_dt_p

  !---------------------------------------------!
  !       Store old particle coordinates        !
  !   (important for bouncing from the walls)   !
  !---------------------------------------------!
  do k = 1, Swarm % n_particles
    Part => Swarm % Particle(k)
    Part % x_o = Part % x_n
    Part % y_o = Part % y_n
    Part % z_o = Part % z_n
  end do

  ! Store old values of smoothed vof function
  do k = 1, Swarm % n_particles
    Part => Swarm % Particle(k)
    Part % smooth_o = Part % smooth_n
  end do

  !----------------------------------!
  !                                  !
  !   Browse through all particles   !
  !                                  !
  !----------------------------------!
  do ss = 1, Swarm % n_sub_steps

    n_parts_in_buffers = 0

    do k = 1, Swarm % n_particles

      ! Take aliases for the particle
      Part      => Swarm % Particle(k)
      deposited => Part  % deposited
      escaped   => Part  % escaped
      trapped   => Part  % trapped

      !--------------------------------------------------!
      !   If particle is neither deposited nor escaped   !
      !--------------------------------------------------!
      if(.not. deposited .and. .not. escaped .and. .not. trapped) then

        ! If particle is in this processor, carry on with it
        if(Part % proc .eq. This_Proc()) then

          ! Compute velocity at the particle, and move it
          ! (also calls Bounce_Particle)
          call Swarm % Move_Particle(k)

          ! Handle its interaction with boundaries
          call Swarm % Bounce_Particle(k)

          ! Handle its interaction with surface
          call Swarm % Trap_Particle(k)

          ! Calling particle forces subroutine to ...
          ! ... compute the forces on each particle and store it
          call Swarm % Particle_Forces(k)

          ! Calling the nearest cell subroutine to find the ...
          ! ... nearest cell for each particle and stores it
          call Part % Find_Nearest_Cell(n_parts_in_buffers)

          ! Calling the nearest node subroutine to find the ...
          ! ... nearest node for each particle and stores it
          call Part % Find_Nearest_Node()

          ! First check if it didn't escape through periodicity
          call Swarm % Check_Periodicity(k, n_parts_in_buffers)

          ! Gathering Swarm statistics
          call Swarm % Calculate_Particles_Mean(k, n_stat_p)

        end if  ! in this processor
      end if    ! deposited or escaped
    end do      ! through particles

    ! Exchange particles for parallel version; if needed
    call Global % Sum_Int(n_parts_in_buffers)
    if(n_parts_in_buffers > 0) then
      call Swarm % Exchange_Particles()
    end if

  end do        ! through sub-steps

  !---------------------------------!
  !   Move particles on a surface   !
  !---------------------------------!
  n_parts_in_buffers = 0

  do k = 1, Swarm % n_particles

    ! Take aliases for the particle
    Part    => Swarm % Particle(k)
    trapped => Part  % trapped

    if(trapped) then

      ! If particle is in this processor, carry on with it
      if(Part % proc .eq. This_Proc()) then
        call Swarm % Move_Trapped(k)

        ! It might have moved to a new cell
        call Part % Find_Nearest_Cell(n_parts_in_buffers)
        call Part % Find_Nearest_Node()

      end if

    end if    ! if trapped

  end do      ! through particles

  ! Exchange particles for parallel version; if needed
  call Global % Sum_Int(n_parts_in_buffers)
  if(n_parts_in_buffers > 0) then
    call Swarm % Exchange_Particles()
  end if

  !-----------------------------------!
  !   Print some data on the screen   !
  !-----------------------------------!
  call Swarm % Print_Swarm_Statistics()

  end subroutine
