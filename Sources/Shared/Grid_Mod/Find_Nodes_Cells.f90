!==============================================================================!
  subroutine Grid_Mod_Find_Nodes_Cells(grid)
!------------------------------------------------------------------------------!
!   Cells around each node are needed for Lagrangian particle tracking.        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, ln, n, run, s1, s2, n1, n2, ln1, ln2, nc1, nc2, nu
  integer              :: max_n_cells
  integer, allocatable :: cell_list(:)
  real                 :: x1, x2, y1, y2, z1, z2

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
  !   and allocate memory in run 1, store them in run 2    !
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
      allocate(cell_list(max_n_cells*2));  cell_list(:) = 0
    end if

  end do

  !-------------------------------------------------------------------!
  !   Connect nodes' cell lists through periodic boundaries as well   !
  !-------------------------------------------------------------------!

  ! Run through three possible periodic directions.
  ! (If you do it less, some cells might get missed
  !  during the browsing and gathering procedure)
  do run = 1, 3
    do s1 = 1, grid % n_faces
      s2 = grid % faces_s(s1)

      ! Take only faces with shadows - ony those have shadow nodes
      if(s2 > 0) then
        do ln1 = 1, grid % faces_n_nodes(s1)
          n1 = grid % faces_n(ln1, s1)
          x1 = grid % xn(n1)
          y1 = grid % yn(n1)
          z1 = grid % zn(n1)
          do ln2 = 1, grid % faces_n_nodes(s2)
            n2 = grid % faces_n(ln2, s2)
            x2 = grid % xn(n2) + grid % dx(s2)
            y2 = grid % yn(n2) + grid % dy(s2)
            z2 = grid % zn(n2) + grid % dz(s2)

            ! Nodes are matching, add cells to each other's list
            if(Math_Mod_Distance_Squared(x1, y1, z1, x2, y2, z2) < PICO) then

              ! Gather list of cells for both nodes and make a unique sore
              nc1 = grid % nodes_n_cells(n1)  ! number of cells arund node 1
              nc2 = grid % nodes_n_cells(n2)  ! number of cells arund node 2
              cell_list(    1:nc1    ) = grid % nodes_c(1:nc1, n1)
              cell_list(nc1+1:nc1+nc2) = grid % nodes_c(1:nc2, n2)
              call Sort_Mod_Unique_Int(cell_list(1:nc1+nc2), nu)

              ! Copy unique list of cells to both nodes' lists
              grid % nodes_n_cells(n1) = nu
              grid % nodes_n_cells(n2) = nu
              grid % nodes_c(1:nu, n1) = cell_list(1:nu)
              grid % nodes_c(1:nu, n2) = cell_list(1:nu)
            end if
          end do
        end do
      end if
    end do
  end do

! ! Check #1, save those nodes' cells
! do c = 1, grid % n_cells - grid % comm % n_buff_cells
!   do ln = 1, grid % cells_n_nodes(c)  ! local node number
!     n = grid % cells_n(ln, c)         ! global node number
!     write(600+this_proc, '(a,i2,a,36i6)')  &
!           'nc=', grid % nodes_n_cells(n),  &
!           ' cell list=', grid % nodes_c(1:max_n_cells, n)
!     write(600+this_proc, '(a,36i6)')  &
!           '           proc=',  &
!           grid % comm % cell_proc(grid % nodes_c(1:max_n_cells, n))
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
