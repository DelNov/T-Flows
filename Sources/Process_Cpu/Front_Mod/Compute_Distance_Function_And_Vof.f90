!==============================================================================!
  subroutine Compute_Distance_Function_And_Vof(Front, cell_dist, vof)
!------------------------------------------------------------------------------!
!>  This subroutine calculates the distance from a surface and stores it
!>  in parameter 'cell_dist'.  It computes a distance function and updates
!>  (corrects) the VOF field for each cell in a computational domain.
!>  At the time of writing this documentation, this subroutine was not used.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Establishing connections to the real cell data.                          !
!   * Iterating through vertices to determine the nearest cells and create a   !
!     list of surrounding nodes and cells.                                     !
!   * Comparing distance measures to find the minimum and maximum distances.   !
!   * Transforming these distances into a specific function and finally        !
!     computing the VOF based on this function.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type), target :: Front      !! parent class
  type(Var_Type)            :: cell_dist  !! variable which stores distances
  type(Var_Type)            :: vof        !! variable representing VOF
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer  :: Grid
  real                      :: dist
  real                      :: xc, yc, zc, xs, ys, zs, min_d, max_d, eps
  integer                   :: nc, nn, nv, n, c, d, v, e
  integer                   :: i_ele, i_nod, j_nod, i_cel
  integer                   :: n_cnt, node_list(512)
  integer                   :: c_cnt, cell_list(512)
  real, contiguous, pointer :: phi_c(:)
!==============================================================================!

  call Work % Connect_Real_Cell(phi_c)

  ! Take some aliases
  Grid => Front % pnt_grid
  nv = Front % n_verts
  nc = Grid % n_cells
  nn = Grid % n_nodes

  ! Initialize distance to some large value
  phi_c(1:nc) = MEGA
  min_d =  MEGA
  max_d = -MEGA

  !---------------------------------!
  !   Browse through all vertices   !
  !---------------------------------!
  do v = 1, nv

    ! Take nearest cell
    d = Front % Vert(v) % cell  ! nearest cell

    ! Make a list of nodes surrounding the nearest cell
    n_cnt = 0  ! initialize node count
    do j_nod = 1, abs(Grid % cells_n_nodes(d))
      n = Grid % cells_n(j_nod, d)

      do i_cel = 1, Grid % nodes_n_cells(n)
        c = Grid % nodes_c(i_cel, n)
        do i_nod = 1, abs(Grid % cells_n_nodes(c))
          n_cnt = n_cnt + 1
          node_list(n_cnt) = Grid % cells_n(i_nod, c)
        end do
      end do
    end do

    ! Compress the list of nodes
    call Sort % Unique_Int(node_list(1:n_cnt), n_cnt)

    ! Make a list of cells surrounding the stored nodes
    c_cnt = 0  ! initialize cell count
    do i_nod = 1, n_cnt
      do i_cel = 1, Grid % nodes_n_cells( node_list(i_nod) )
        c_cnt = c_cnt + 1
        cell_list(c_cnt) = Grid % nodes_c(i_cel, node_list(i_nod))
      end do
    end do
    call Sort % Unique_Int(cell_list(1:c_cnt), c_cnt)

    ! Match first and second neighbour cells with all elements around the node
    do i_ele = 1, Front % Vert(v) % nne
      e  = Front % Vert(v) % e(i_ele)

      ! Center of the sphere
      xs = Front % Elem(e) % xc
      ys = Front % Elem(e) % yc
      zs = Front % Elem(e) % zc

      do i_cel = 1, c_cnt
        c = cell_list(i_cel)

        xc = Grid % xc(c)
        yc = Grid % yc(c)
        zc = Grid % zc(c)

        dist = Math % Distance(xc, yc, zc, xs, ys, zs)  &
             - 1.0 / Front % Elem(e) % curv

        if(abs(dist) < abs(phi_c(c))) then
          phi_c(c) = dist
          min_d = min(phi_c(c), min_d)
          max_d = max(phi_c(c), max_d)
        end if
      end do
    end do

  end do  ! through vertices

  call Global % Min_Real(min_d)
  call Global % Max_Real(max_d)

  !----------------------------------------------------------------------!
  !   Transform the distances into what Mijail calls \phi in his paper   !
  !----------------------------------------------------------------------!
  do c = 1, nc
    cell_dist % n(c) = phi_c(c)
    if( Math % Approx_Real(phi_c(c), MEGA) ) then
      if(vof % n(c) < 0.5) then
        cell_dist % n(c) = min_d
      else
        cell_dist % n(c) = max_d
      end if
    end if
  end do

  !---------------------------------------------------!
  !   Finally, calculate VOF from distance function   !
  !---------------------------------------------------!
  eps = sqrt(2.0) * 1.0 / 64  ! hard-coded epsilon
  do c = 1, nc
    if(cell_dist % n(c) < -eps) then
      vof % n(c) = 0.0
    else if(cell_dist % n(c) > eps) then
      vof % n(c) = 1.0
    else
      vof % n(c) = 0.5 * (  1.0                                          &
                          + cell_dist % n(c) / eps                       &
                          + 1.0/PI * sin(PI * cell_dist % n(c) / eps) )
    end if
  end do

  call Work % Disconnect_Real_Cell(phi_c)

  end subroutine
