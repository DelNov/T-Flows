!==============================================================================!
  subroutine Swarm_Mod_Find_Nearest_Node(flow, swarm, k)
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
  type(Field_Type), target :: flow
  type(Swarm_Type), target :: swarm
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: grid
  type(Particle_Type), pointer :: part
  integer                      :: ln, n        ! local node and node
  integer                      :: cn           ! closest node
  integer                      :: k            ! particle count 
  integer                      :: n_part       ! number of particles
  real                         :: xn, yn, zn   ! Node coordinates
  real                         :: xc, yc, zc   ! cell center coordinates
  real                         :: dn_sq        ! distance squared
  real                         :: min_dn       ! minimum distance computed
  real,                pointer :: xp, yp, zp   ! coordinates of particle
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  part => swarm % particles(k)

  !-------------------------------------------------!
  !   Browse through all cells as another example   !
  !-------------------------------------------------!
  min_dn = HUGE

  ! Browse through nodes of particle's cell
  do ln = 1, grid % cells_n_nodes(part % cell)  ! local node number
    n = grid % cells_n(ln, part % cell)         ! global node number

    ! Take node coordinate
    xn = grid % xn(n)
    yn = grid % yn(n)
    zn = grid % zn(n)

    ! Distance squared from the particle to cell's node
    dn_sq = (xn - part % x)**2  &
          + (yn - part % y)**2  &
          + (zn - part % z)**2

    ! Finding the closest node of those 4 of the closest cell
    if(dn_sq < min_dn) then
      min_dn = dn_sq  ! new minimum distance
      cn = n          ! cn is the closest node index
    end if

  end do

  part % node = cn

  end subroutine
