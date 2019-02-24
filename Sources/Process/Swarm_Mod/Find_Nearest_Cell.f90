!==============================================================================!
  subroutine Swarm_Mod_Find_Nearest_Cell(swarm, k)
!------------------------------------------------------------------------------!
!   Finds a cell closest to a particle.                                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
  integer                  :: k      ! particle number
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: grid
  type(Particle_Type), pointer :: part
  integer                      :: c, cc       ! cell and closest cell
  integer                      :: n_part      ! particle number
  real                         :: xc, yc, zc  ! cell center coordinates
  real                         :: d_sq        ! distance squared
  real                         :: min_d       ! minimum distance computed
  logical,             pointer :: reflected   ! reflection BC type
!==============================================================================!

  ! Take aliases
  grid      => swarm % pnt_grid
  part      => swarm % particle(k)
  reflected => part  % reflected

  !-------------------------------------------------!
  !                                                 !
  !   Browse through all cells as another example   !
  !                                                 !
  !-------------------------------------------------!
  min_d = HUGE

  !-------------!
  !   Trap BC   !
  !-------------!
  if (.not. reflected) then
    if(.not. part % escaped .and. .not. part % deposited) then 

      do c = 1, grid % n_cells

        ! Take cell centre
        xc = grid % xc(c)
        yc = grid % yc(c)
        zc = grid % zc(c)

      ! Distance squared from the particle to cell centre
        d_sq = (xc - part % x)**2  &
             + (yc - part % y)**2  &
             + (zc - part % z)**2

        if(d_sq < min_d) then 
          min_d = d_sq  ! new minimum distance
          cc = c        ! cc is the closest cell index
        end if

      end do

      ! save the last value of cc (the index of closest cell) 
      part % cell = cc

      ! Calling the nearest node subroutine to find the ...
      ! ... nearest node for each particle and stores it
      call Swarm_Mod_Find_Nearest_Node(swarm, k)

    else
      reflected = .true.  ! end because the particle either escaped or deposited
    end if
  end if

  !------------------!
  !   Reflected BC   !
  !------------------!
  if (reflected) then
    if(.not. part % escaped) then

      do c = 1, grid % n_cells

        ! Take cell centre
        xc = grid % xc(c)
        yc = grid % yc(c)
        zc = grid % zc(c)

        ! Distance squared from the particle to cell centre
        d_sq = (xc - part % x)**2  &
             + (yc - part % y)**2  &
             + (zc - part % z)**2

        if(d_sq < min_d) then
          min_d = d_sq  ! new minimum distance
          cc = c        ! cc is the closest cell index
        end if
      end do

      ! saves the last value of cc (the index of closest cell) 
      part % cell = cc

      ! Calling the nearest node subroutine to find the ...
      ! ... nearest node for each particle and stores it
      call Swarm_Mod_Find_Nearest_Node(swarm, k)

    else
      reflected = .false.
    end if
  end if

  end subroutine
