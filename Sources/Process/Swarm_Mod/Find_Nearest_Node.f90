!==============================================================================!
  subroutine Swarm_Mod_Find_Nearest_Node(swarm, k)
!------------------------------------------------------------------------------!
!   Finds a node closest to a particle.                                        !
!   Important: it assumes that the closest cell has been found!                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
  integer                  :: k      ! particle number
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),    pointer :: flow
  type(Grid_Type),     pointer :: grid
  type(Particle_Type), pointer :: part
  integer                      :: ln, n        ! local node and node
  integer                      :: cn           ! closest node
  real                         :: xn, yn, zn   ! Node coordinates
  real                         :: dn_sq        ! distance squared
  real                         :: min_dn       ! minimum distance computed
!==============================================================================!

  ! Take aliases
  flow => swarm % pnt_flow
  grid => swarm % pnt_grid
  part => swarm % particle(k)

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
    dn_sq = (xn - part % x_n)**2  &
          + (yn - part % y_n)**2  &
          + (zn - part % z_n)**2

    ! Finding the closest node of those 4 of the closest cell
    if(dn_sq < min_dn) then
      min_dn = dn_sq  ! new minimum distance
      cn = n          ! cn is the closest node index
    end if
  end do

  part % node = cn

  end subroutine
