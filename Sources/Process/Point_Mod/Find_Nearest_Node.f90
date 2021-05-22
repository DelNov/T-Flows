!==============================================================================!
  subroutine Find_Nearest_Node(Point)
!------------------------------------------------------------------------------!
!   Finds a node closest to the point                                          !
!   Important: it assumes that the closest cell has been found!                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Point_Type), target :: Point
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: i_nod, n     ! local node and node
  integer                  :: cn           ! closest node
  real                     :: xn, yn, zn   ! Node coordinates
  real                     :: dn_sq        ! distance squared
  real                     :: min_dn       ! minimum distance computed
!==============================================================================!

  ! Take aliases
  Grid => Point % pnt_grid

  !-------------------------------------------------!
  !   Browse through all cells as another example   !
  !-------------------------------------------------!
  min_dn = HUGE

  ! Browse through nodes of particle's cell
  do i_nod = 1, abs(Grid % cells_n_nodes(Point % cell))  ! local node number
    n = Grid % cells_n(i_nod, Point % cell)              ! global node number

    ! Take node coordinate
    xn = Grid % xn(n)
    yn = Grid % yn(n)
    zn = Grid % zn(n)

    ! Distance squared from the particle to cell's node
    dn_sq = (xn - Point % x_n)**2  &
          + (yn - Point % y_n)**2  &
          + (zn - Point % z_n)**2

    ! Finding the closest node of those 4 of the closest cell
    if(dn_sq < min_dn) then
      min_dn = dn_sq  ! new minimum distance
      cn = n          ! cn is the closest node index
    end if
  end do

  Point % node = cn

  end subroutine
