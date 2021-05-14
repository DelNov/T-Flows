!==============================================================================!
  subroutine Find_Nearest_Node(Front, v)
!------------------------------------------------------------------------------!
!   Finds a node closest to a vertex.                                          !
!   Important: it assumes that the closest cell has been found!                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type), target :: Front
  integer                   :: v      ! vertex number
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  type(Vert_Type),  pointer :: vert
  integer                   :: i_nod, n     ! local node and node
  integer                   :: cn           ! closest node
  real                      :: xn, yn, zn   ! node coordinates
  real                      :: dn_sq        ! distance squared
  real                      :: min_dn       ! minimum distance computed
!==============================================================================!

  ! Take aliases
  flow => Front % pnt_flow
  grid => Front % pnt_grid
  vert => Front % vert(v)

  !-------------------------------------------------!
  !   Browse through all cells as another example   !
  !-------------------------------------------------!
  min_dn = HUGE

  ! Browse through nodes of vertex's cell
  do i_nod = 1, grid % cells_n_nodes(vert % cell)  ! local node number
    n = grid % cells_n(i_nod, vert % cell)         ! global node number

    ! Take node coordinate
    xn = grid % xn(n)
    yn = grid % yn(n)
    zn = grid % zn(n)

    ! Distance squared from the vertex to cell's node
    dn_sq = (xn - vert % x_n)**2  &
          + (yn - vert % y_n)**2  &
          + (zn - vert % z_n)**2

    ! Finding the closest node of those 4 of the closest cell
    if(dn_sq < min_dn) then
      min_dn = dn_sq  ! new minimum distance
      cn = n          ! cn is the closest node index
    end if

  end do

  vert % node = cn

  end subroutine
