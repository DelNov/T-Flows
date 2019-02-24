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
  integer                      :: c, lc, cc   ! cell, local and closest cell
  real                         :: xc, yc, zc  ! cell center coordinates
  real                         :: dc_sq       ! distance squared
  real                         :: min_dc      ! minimum distance computed
!==============================================================================!

  ! Take aliases
  grid => swarm % pnt_grid
  part => swarm % particle(k)

  ! Initialize minimum distance
  min_dc = HUGE

  !-----------------------------------------------------------!
  !   Closest node is known, browse through cells around it   !
  !-----------------------------------------------------------!
  if(part % node .ne. 0) then

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
        cc = c          ! cc is the closest cell index
      end if

    end do

  !--------------------------------------------------------!
  !   Closest node is not known, brows through all cells   !
  !--------------------------------------------------------!
  else

    do c = 1, grid % n_cells

      ! Take cell centre
      xc = grid % xc(c)
      yc = grid % yc(c)
      zc = grid % zc(c)

      ! Distance squared from the particle to cell centre
      dc_sq = (xc - part % x)**2  &
            + (yc - part % y)**2  &
            + (zc - part % z)**2

      if(dc_sq < min_dc) then
        min_dc = dc_sq  ! new minimum distance
        cc = c          ! cc is the closest cell index
      end if

    end do

  end if

  ! Save the last value of cc (the index of closest cell) 
  part % cell = cc

  end subroutine
