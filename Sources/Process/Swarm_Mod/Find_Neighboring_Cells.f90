!==============================================================================!
  subroutine Swarm_Mod_Find_Neighboring_Cells(flow, swarm, k)
!------------------------------------------------------------------------------!
!   Finds a cell closest to a particle.                                        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod, only: HUGE
  use Field_Mod, only: Field_Type
  use Grid_Mod,  only: Grid_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: flow
  type(Swarm_Type),    target    :: swarm
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: grid
  type(Particle_Type), pointer :: part
  integer                      :: lc, c        ! local cell and cell
  integer                      :: n            ! node counter
  integer                      :: cn           ! closest node
  integer                      :: cc           ! closest sell
  integer                      :: k            ! particle count 
  integer                      :: n_part       ! number of particles
  real                         :: xn, yn, zn   ! Node coordinates
  real                         :: xc, yc, zc   ! cell center coordinates
  real                         :: dc_sq        ! distance squared
  real                         :: min_dc       ! minimum distance computed
  logical,             pointer :: escaped
  logical,             pointer :: deposited
  logical,             pointer :: trapped      ! trap       BC type
  logical,             pointer :: reflected    ! reflection BC type
!==============================================================================!

  ! Take aliases
  grid         => flow % pnt_grid
  part         => swarm % particles(k)
  escaped      => part  % escaped
  deposited    => part  % deposited
  reflected    => part % reflected
  trapped      => part % trapped

  !-------------------------------------------------!
  !   Browse through all cells as another example   !
  !-------------------------------------------------!
  min_dc = HUGE

  !----------------!
  !   Trapped BC   !
  !----------------!
  if (trapped) then
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
      call Swarm_Mod_Find_Nearest_Node(flow, swarm, k)

      ! Calling interpolate velocity subroutine ...
      ! ... to update velocity of each particle
      call Swarm_Mod_Interpolate_Velocity(flow, swarm, k)

      ! Calling particle forces subroutine to ...
      ! ... compute the forces on each particle and store it
      call Swarm_Mod_Particle_Forces(flow, swarm, k)

    else
      trapped = .false. ! End because the particle either escaped or deposited
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
      call Swarm_Mod_Find_Nearest_Node(flow, swarm, k)

      ! Calling interpolate velocity subroutine ...
      ! ... to update velocity of each particle
      call Swarm_Mod_Interpolate_Velocity(flow, swarm, k)

      ! Calling particle forces subroutine to ...
      ! ... compute the forces on each particle and store it
      call Swarm_Mod_Particle_Forces(flow, swarm, k)

    else
      reflected = .false.
    end if

  end if

  end subroutine
