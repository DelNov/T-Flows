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
  integer                      :: lc, c        ! local cell and cell
  integer                      :: n            ! node counter
  integer                      :: cn           ! closest node
  integer                      :: cc           ! closest sell
  integer                      :: n_part       ! number of particles
  real                         :: xn, yn, zn   ! Node coordinates
  real                         :: xc, yc, zc   ! cell center coordinates
  real                         :: dc_sq        ! distance squared
  real                         :: min_dc       ! minimum distance computed
  logical,             pointer :: escaped
  logical,             pointer :: deposited
  logical,             pointer :: reflected    ! reflection BC type
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
    reflected => part  % reflected

    !-------------------------------------------------!
    !   Browse through all cells as another example   !
    !-------------------------------------------------!
    min_dc = HUGE

    !----------------!
    !   Trapped BC   !
    !----------------!
    if (.not. reflected) then
      if(.not. escaped .and. .not. deposited) then 

        !------------------------------------------------!
        !   Now we're going to use that node index ...   !
        !   ... to search for the 8 neighboring cells    !
        !------------------------------------------------!
        do n = 1, grid % n_nodes
          do lc = 1, grid % nodes_n_cells(part % node)   ! local cell number
            c = grid % nodes_c(lc, part % node)          ! global cell number

            xc = grid % xc(c)
            yc = grid % yc(c)
            zc = grid % zc(c)

            ! Distance squared from the particle to cell centre
            dc_sq = (xc - part % x)**2  &
                  + (yc - part % y)**2  &
                  + (zc - part % z)**2

            ! Finding the closest cell of those 8 neighbors
            if(dc_sq < min_dc) then
              min_dc = dc_sq  ! new minimum distance
              cc = c        ! cc is the closest cell index
            end if

          end do
        end do

        ! Closest cell of the 8 neighboring to the node
        part % cell = cc

        ! Calling the nearest node subroutine to find the ...
        ! ... nearest node for each particle and stores it
        call Swarm_Mod_Find_Nearest_Node(swarm, k)

        ! Calling interpolate velocity subroutine ...
        ! ... to update velocity of each particle
        call Swarm_Mod_Interpolate_Velocity(swarm, k)

        ! Calling particle forces subroutine to ...
        ! ... compute the forces on each particle and store it
        call Swarm_Mod_Particle_Forces(swarm, k)

      else
        reflected = .true. ! End because the particle either escaped or deposited
      end if
    end if

    !-------------------!
    !   Reflection BC   !
    !-------------------!
    if (reflected) then
      if(.not. escaped) then 

        !------------------------------------------------!
        !   Now we're going to use that node index ...   !
        !   ... to search for the 8 neighboring cells    !
        !------------------------------------------------!
        do n = 1, grid % n_nodes
          do lc = 1, grid % nodes_n_cells(part % node)    ! local cell number
             c = grid % nodes_c(lc, part % node)          ! global cell number

             xc = grid % xc(c)
             yc = grid % yc(c)
             zc = grid % zc(c)

            ! Distance squared from the particle to cell centre
            dc_sq = (xc - part % x)**2  &
                  + (yc - part % y)**2  &
                  + (zc - part % z)**2

            ! Finding the closest cell of those 8 neighbors
             if(dc_sq < min_dc) then
               min_dc = dc_sq  ! new minimum distance
               cc = c        ! cc is the closest cell index
             end if

          end do
        end do

        ! Closest cell of the 8 neighboring to the node
        part % cell = cc

        ! Calling the nearest node subroutine to find the ...
        ! ... nearest node for each particle and stores it
        call Swarm_Mod_Find_Nearest_Node(swarm, k)

        ! Calling interpolate velocity subroutine ...
        ! ... to update velocity of each particle
        call Swarm_Mod_Interpolate_Velocity(swarm, k)

        ! Calling particle forces subroutine to ...
        ! ... compute the forces on each particle and store it
        call Swarm_Mod_Particle_Forces(swarm, k)

      else
        reflected = .false.
      end if

    end if

  end do

  end subroutine
