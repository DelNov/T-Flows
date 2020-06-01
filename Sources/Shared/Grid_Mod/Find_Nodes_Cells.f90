!==============================================================================!
  subroutine Grid_Mod_Find_Nodes_Cells(grid)
!------------------------------------------------------------------------------!
!   Cells around each node are needed for Lagrangian particle tracking.        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, ln, n, run   ! counters
  integer :: max_n_cells     ! max number of cells surrounding a node

! Just for checking, erase this later
! integer :: lc
! real    :: xn, xc, yn, yc, zn, zc, max_del, min_del, del
!==============================================================================!

  ! Allocate memory for node coordinates
  n = grid % n_nodes
  allocate(grid % nodes_n_cells(1:n))
  grid % nodes_n_cells(:) = 0

  !--------------------------------------------------------!
  !   Find maximum number of cells surrounding each node   !
  !   (That includes sells inside and on the boundaries)   !
  !--------------------------------------------------------!

  do run = 1, 2
    grid % nodes_n_cells(:) = 0  ! re-initialize the cell count

    do c = -grid % n_bnd_cells, grid % n_cells
      do ln = 1, grid % cells_n_nodes(c)  ! local node number
        n = grid % cells_n(ln, c)         ! global node number

        ! Increase number of cells surrounding the this node by one
        grid % nodes_n_cells(n) = grid % nodes_n_cells(n) + 1

        ! Store the current cell in the second run
        if(run .eq. 2) then
          grid % nodes_c(grid % nodes_n_cells(n), n) = c
        end if
      end do
    end do

    ! Allocate memory for cells surrounding each node
    if(run .eq. 1) then
      max_n_cells = maxval(grid % nodes_n_cells)
      allocate(grid % nodes_c(1:max_n_cells, 1:grid % n_nodes))
    end if

  end do

! ! Check #1, save those nodes' cells
! do c = 1, grid % n_cells - grid % comm % n_buff_cells
!   do ln = 1, grid % cells_n_nodes(c)  ! local node number
!     n = grid % cells_n(ln, c)         ! global node number
!     write(500+this_proc, '(36i6)') grid % nodes_n_cells(n)
!     write(600+this_proc, '(36i6)') grid % nodes_c(1:max_n_cells, n)
!   end do
! end do

! ! Check #2, erase this later
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
