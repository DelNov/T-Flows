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
!==============================================================================!

  ! Take aliases
  grid => swarm % pnt_grid
  part => swarm % particle(k)

  ! Initialize minimum distance
  min_d = HUGE

  !-------------------------------------------------!
  !   Browse through all cells as another example   !
  !-------------------------------------------------!
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

  ! Save the last value of cc (the index of closest cell) 
  part % cell = cc

  end subroutine
