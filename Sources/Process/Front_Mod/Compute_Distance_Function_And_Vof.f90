!==============================================================================!
  subroutine Front_Mod_Compute_Distance_Function_And_Vof(front, cell_dist, vof)
!------------------------------------------------------------------------------!
!   Computes distance from a surface                                           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: phi_c  => r_cell_01  ! cell values of phi
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Front_Type), target :: front
  type(Var_Type)           :: cell_dist
  type(Var_Type)           :: vof
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: grid
  real                     :: dist
  real                     :: xc, yc, zc, xs, ys, zs, min_d, max_d, eps
  integer                  :: nc, nn, nv, n, c, d, v, e
  integer                  :: i_ele, i_nod, j_nod, i_cel
  integer                  :: n_cnt, node_list(512)
  integer                  :: c_cnt, cell_list(512)
!==============================================================================!

  ! Take some aliases
  grid => front % pnt_grid
  nv = front % n_verts
  nc = grid % n_cells
  nn = grid % n_nodes

  ! Initialize distance to some large value
  phi_c(1:nc) = MEGA
  min_d =  MEGA
  max_d = -MEGA

  !---------------------------------!
  !   Browse through all vertices   !
  !---------------------------------!
  do v = 1, nv

    ! Take nearest cell
    d = front % vert(v) % cell  ! nearest cell

    ! Make a list of nodes surrounding the nearest cell
    n_cnt = 0  ! initialize node count
    do j_nod = 1, grid % cells_n_nodes(d)
      n = grid % cells_n(j_nod, d)

      do i_cel = 1, grid % nodes_n_cells(n)
        c = grid % nodes_c(i_cel, n)
        do i_nod = 1, grid % cells_n_nodes(c)
          n_cnt = n_cnt + 1
          node_list(n_cnt) = grid % cells_n(i_nod, c)
        end do
      end do
    end do

    ! Compress the list of nodes
    call Sort % Unique_Int(node_list(1:n_cnt), n_cnt)

    ! Make a list of cells surrounding the stored nodes
    c_cnt = 0  ! initialize cell count
    do i_nod = 1, n_cnt
      do i_cel = 1, grid % nodes_n_cells( node_list(i_nod) )
        c_cnt = c_cnt + 1
        cell_list(c_cnt) = grid % nodes_c(i_cel, node_list(i_nod))
      end do
    end do
    call Sort % Unique_Int(cell_list(1:c_cnt), c_cnt)

    ! Match first and second neighbour cells with all elements around the node
    do i_ele = 1, front % vert(v) % nne
      e  = front % vert(v) % e(i_ele)

      ! Center of the sphere
      xs = front % elem(e) % xc
      ys = front % elem(e) % yc
      zs = front % elem(e) % zc

      do i_cel = 1, c_cnt
        c = cell_list(i_cel)

        xc = grid % xc(c)
        yc = grid % yc(c)
        zc = grid % zc(c)

        dist = Math_Mod_Distance(xc, yc, zc, xs, ys, zs)  &
                                 - 1.0 / front % elem(e) % curv

        if(abs(dist) < abs(phi_c(c))) then
          phi_c(c) = dist
          min_d = min(phi_c(c), min_d)
          max_d = max(phi_c(c), max_d)
        end if
      end do
    end do

  end do  ! through vertices

  call Comm_Mod_Global_Min_Real(min_d)
  call Comm_Mod_Global_Max_Real(max_d)

  !----------------------------------------------------------------------!
  !   Transform the distances into what Mijail calls \phi in his paper   !
  !----------------------------------------------------------------------!
  do c = 1, nc
    cell_dist % n(c) = phi_c(c)
    if( Math_Mod_Approx_Real(phi_c(c), MEGA) ) then
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

  end subroutine
