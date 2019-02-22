!==============================================================================!
  subroutine Grid_Mod_Find_Nodes_Cells(grid)
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod, only: HUGE
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: lc, c, ln, n   ! locel cell, cell, local node and node
  integer :: max_n_cells    ! max number of cells surrounding a node

! Just for checking, erase this later
  real :: xn, xc, yn, yc, zn, zc, max_del, min_del, del
!==============================================================================!

  ! Allocate memory for node coordinates
  n = grid % n_nodes
  allocate(grid % nodes_n_cells(1:n))
  grid % nodes_n_cells(:) = 0

  !--------------------------------------------------------!
  !   Find maximum number of cells surrounding each node   !
  !     (use only inside cells at this point in time)      !
  !--------------------------------------------------------!
  do c = 1, grid % n_cells
    do ln = 1, grid % cells_n_nodes(c)  ! local node number
      n = grid % cells_n(ln, c)         ! global node number

      ! Increase number of cells surrounding the this node by one
      grid % nodes_n_cells(n) = grid % nodes_n_cells(n) + 1
    end do
  end do

  max_n_cells = maxval(grid % cells_n_nodes)

  ! Allocate memory for cells surrounding each node
  n = grid % n_nodes
  allocate(grid % nodes_c(1:max_n_cells, 1:n))

  !----------------------------------------------------------!
  !   Now you can really store the cells surrounding nodes   !
  !----------------------------------------------------------!
  grid % nodes_n_cells(:) = 0  ! re-initialize the cell count

  do c = 1, grid % n_cells
    do ln = 1, grid % cells_n_nodes(c)  ! local node number
      n = grid % cells_n(ln, c)         ! global node number

      ! Increase number of cells surrounding the this node by one ...
      grid % nodes_n_cells(n) = grid % nodes_n_cells(n) + 1

      ! ... and store the current cell
      grid % nodes_c(grid % nodes_n_cells(n), n) = c
    end do
  end do

! ! Just for checking, erase this later
! max_del = -HUGE
! min_del = +HUGE
!
! do n = 1, grid % n_nodes
!   do lc = 1, grid % nodes_n_cells(n)  ! local cell number
!     c = grid % nodes_c(lc, n)         ! global cell number
!
!     xn = grid % xn(n)
!     yn = grid % yn(n)
!     zn = grid % zn(n)
!
!     xc = grid % xc(c)
!     yc = grid % yc(c)
!     zc = grid % zc(c)
!
!     del = sqrt( (xn-xc)**2 + (yn-yc)**2 + (zn-zc)**2 )
!
!     min_del = min(del, min_del)
!     max_del = max(del, max_del)
!   end do
! end do
!
! print *, '# Checking: min and max del: ', min_del, max_del
! stop

  end subroutine
