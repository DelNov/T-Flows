!==============================================================================!
  subroutine Find_Nodes_Cells(Grid)
!------------------------------------------------------------------------------!
!>  This subroutine identifies and lists all cells surrounding each node
!>  within the computational grid. It is particularly important for Lagrangian
!>  particle tracking, as it provides the necessary spatial information to
!>  accurately trace particle paths through the grid.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initialization:                                                          !
!     - Allocates memory for node-to-cell mapping and initializes the count    !
!       of surrounding cells for each node.                                    !
!   * Cell Mapping:                                                            !
!     - Iteratively maps each cell to its nodes, keeping track of the number   !
!       of cells associated with each node.                                    !
!     - In a two-run approach, first determines the number of cells per node   !
!       and then stores the actual cells in the second run.                    !
!   * Periodic Boundary Consideration:                                         !
!     - Extends the mapping to include cells across periodic boundaries,       !
!       ensuring comprehensive coverage of the grid.                           !
!   * Finalization:                                                            !
!     - Concludes the mapping process, updating the nodes-to-cells association !
!       to include boundary faces where needed.                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid  !! computational grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: i_nod, n, run, n1, n2, i_nod1, i_nod2, nc1, nc2, nu
  integer              :: c, c2, s, s1, s2, max_n_cells
  integer, allocatable :: cell_list(:)
  real                 :: x1, x2, y1, y2, z1, z2
! Just for checking, see sections at the end of the procedure
! integer              :: i_cel
! real                 :: xn, xc, yn, yc, zn, zc, max_del, min_del, del
! real, allocatable    :: n_cells(:)
!==============================================================================!

  ! Allocate memory for node coordinates
  n = Grid % n_nodes
  allocate(Grid % nodes_n_cells(1:n))
  Grid % nodes_n_cells(:) = 0

  !--------------------------------------------------------!
  !   Find maximum number of cells surrounding each node   !
  !   and allocate memory in run 1, store them in run 2    !
  !--------------------------------------------------------!

  do run = 1, 2
    Grid % nodes_n_cells(:) = 0  ! re-initialize the cell count

    do c = -Grid % n_bnd_cells, Grid % n_cells
      do i_nod = 1, abs(Grid % cells_n_nodes(c))  ! local node number
        n = Grid % cells_n(i_nod, c)              ! global node number

        ! Increase number of cells surrounding the this node by one
        Grid % nodes_n_cells(n) = Grid % nodes_n_cells(n) + 1

        ! Store the current cell in the second run
        if(run .eq. 2) then
          Grid % nodes_c(Grid % nodes_n_cells(n), n) = c
        end if
      end do
    end do

    ! Allocate memory for cells surrounding each node
    if(run .eq. 1) then
      max_n_cells = maxval(Grid % nodes_n_cells)
      call Global % Max_Int(max_n_cells)
      allocate(Grid % nodes_c(1:max_n_cells, 1:Grid % n_nodes))
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
    do s1 = 1, Grid % n_faces
      s2 = Grid % faces_s(s1)

      ! Take only faces with shadows - ony those have shadow nodes
      if(s2 > 0) then
        do i_nod1 = 1, Grid % faces_n_nodes(s1)
          n1 = Grid % faces_n(i_nod1, s1)
          x1 = Grid % xn(n1)
          y1 = Grid % yn(n1)
          z1 = Grid % zn(n1)
          do i_nod2 = 1, Grid % faces_n_nodes(s2)
            n2 = Grid % faces_n(i_nod2, s2)
            x2 = Grid % xn(n2) + Grid % dx(s2)
            y2 = Grid % yn(n2) + Grid % dy(s2)
            z2 = Grid % zn(n2) + Grid % dz(s2)

            ! Nodes are matching, add cells to each other's list
            if(Math % Distance_Squared(x1, y1, z1, x2, y2, z2) < PICO) then

              ! Gather list of cells for both nodes and make a unique sore
              nc1 = Grid % nodes_n_cells(n1)  ! number of cells arund node 1
              nc2 = Grid % nodes_n_cells(n2)  ! number of cells arund node 2
              cell_list(    1:nc1    ) = Grid % nodes_c(1:nc1, n1)
              cell_list(nc1+1:nc1+nc2) = Grid % nodes_c(1:nc2, n2)
              call Sort % Unique_Int(cell_list(1:nc1+nc2), nu)

              ! Copy unique list of cells to both nodes' lists
              Grid % nodes_n_cells(n1) = nu
              Grid % nodes_n_cells(n2) = nu
              Grid % nodes_c(1:nu, n1) = cell_list(1:nu)
              Grid % nodes_c(1:nu, n2) = cell_list(1:nu)
            end if
          end do
        end do
      end if
    end do
  end do

  !----------------------------------------------!
  !   Store boundary cells to face connections   !
  !----------------------------------------------!
  do s = 1, Grid % n_faces
    c2 = Grid % faces_c(2,s)
    if(c2 < 0) then
      Grid % cells_bnd_face(c2) = s
    end if
  end do

! ! Check #1, save those nodes' cells
! do c = 1, Grid % n_cells - Grid % Comm % n_buff_cells
!   do i_nod = 1, abs(Grid % cells_n_nodes(c))  ! local node number
!     n = Grid % cells_n(i_nod, c)              ! global node number
!     write(600+this_proc, '(a,i2,a,36i6)')  &
!           'nc=', Grid % nodes_n_cells(n),  &
!           ' cell list=', Grid % nodes_c(1:max_n_cells, n)
!     write(600+this_proc, '(a,36i6)')  &
!           '           proc=',  &
!           Grid % Comm % cell_proc(Grid % nodes_c(1:max_n_cells, n))
!   end do
! end do

! ! Check #2
! max_del = -HUGE
! min_del = +HUGE
!
! do n = 1, Grid % n_nodes
!   do i_cel = 1, Grid % nodes_n_cells(n)  ! local cell number
!     c = Grid % nodes_c(i_cel, n)         ! global cell number
!
!     xn = Grid % xn(n)
!     yn = Grid % yn(n)
!     zn = Grid % zn(n)
!
!     xc = Grid % xc(c)
!     yc = Grid % yc(c)
!     zc = Grid % zc(c)
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

! ! Check #3
! allocate(n_cells(Grid % n_nodes));
! n_cells(1:Grid % n_nodes) = Grid % nodes_n_cells(1:Grid % n_nodes)
! call Grid % Save_Debug_Vtu("n-cells",                &
!                            scalar_node = n_cells,    &
!                            scalar_name = "n_cells")
! call Comm_Mod_End
! stop

  end subroutine
