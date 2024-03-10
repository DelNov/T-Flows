!==============================================================================!
  subroutine Find_Nearest_Node(Point)
!------------------------------------------------------------------------------!
!>  This subroutine identifies the node closest to a given point within the
!>  computational grid. It assumes that the closest cell to the point has been
!>  previously determined.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initialization: Prepares variables for distance computation.             !
!   * Node proximity analysis: Calculates the distance to each node within     !
!     the closest cell to identify the nearest node to the point.              !
!   * Closest node determination: Establishes the nearest node based on        !
!     the minimum computed distance.                                           !
!   * Point update: Updates the point's properties with the identified         !
!     closest node.                                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Point_Type), target :: Point  !! point object
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
