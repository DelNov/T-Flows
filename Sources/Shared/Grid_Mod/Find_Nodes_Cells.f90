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
  integer :: lc, c, c2, ln, n, s   ! counters
  integer :: max_n_cells           ! max number of cells surrounding a node

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

  ! Inside cells
  do c = 1, grid % n_cells
    do ln = 1, grid % cells_n_nodes(c)  ! local node number
      n = grid % cells_n(ln, c)         ! global node number

      ! Increase number of cells surrounding the this node by one
      grid % nodes_n_cells(n) = grid % nodes_n_cells(n) + 1
    end do
  end do

  ! Boundary cells
  do s = 1, grid % n_bnd_cells
    c2 = grid % faces_c(2, s)
    if(c2 >= 0) then
      print *, 'PANIC!  Something is very wrong in Find_Nodes_Cells'
    end if
    do ln = 1, grid % faces_n_nodes(s)  ! local face number
      n = grid % faces_n(ln, s)         ! global node number

      ! Increase number of cells surrounding the this node by one
      grid % nodes_n_cells(n) = grid % nodes_n_cells(n) + 1
    end do
  end do

  max_n_cells = maxval(grid % nodes_n_cells)

  ! Allocate memory for cells surrounding each node
  n = grid % n_nodes
  allocate(grid % nodes_c(1:max_n_cells, 1:n))

  !----------------------------------------------------------!
  !   Now you can really store the cells surrounding nodes   !
  !----------------------------------------------------------!
  grid % nodes_n_cells(:) = 0  ! re-initialize the cell count

  ! Inside cells
  do c = 1, grid % n_cells
    do ln = 1, grid % cells_n_nodes(c)  ! local node number
      n = grid % cells_n(ln, c)         ! global node number

      ! Increase number of cells surrounding the this node by one ...
      grid % nodes_n_cells(n) = grid % nodes_n_cells(n) + 1

      ! ... and store the current cell
      grid % nodes_c(grid % nodes_n_cells(n), n) = c
    end do
  end do

  ! Boundary cells
  do s = 1, grid % n_bnd_cells
    c2 = grid % faces_c(2,s)
    do ln = 1, grid % faces_n_nodes(s)  ! local face number
      n = grid % faces_n(ln, s)         ! global node number

      ! Increase number of cells surrounding the this node by one ...
      grid % nodes_n_cells(n) = grid % nodes_n_cells(n) + 1

      ! ... and store the current cell
      grid % nodes_c(grid % nodes_n_cells(n), n) = c2

      ! Also store boundary face for boundary cell
      grid % cells_bnd_face(c2) = s
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
