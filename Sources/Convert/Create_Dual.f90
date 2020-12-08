!==============================================================================!
  subroutine Create_Dual(prim, dual)
!------------------------------------------------------------------------------!
!   Creates a dual mesh from an existing                                       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
  use Save_Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: prim, dual
!----------------------------------[Calling]-----------------------------------!
  integer :: N_Sharp_Corners
  integer :: N_Edges_On_Bnd_Color
  integer :: N_Bnd_Cells_In_Color
  integer :: N_Sharp_Edges
  integer :: N_Nodes_In_Bnd_Color
!-----------------------------------[Locals]-----------------------------------!
  integer              :: bc, e, c, c1, c2, s, n, n1, n2
  integer              :: i, i_nod, j_nod, i_edg
  integer              :: n_bc, d_nn  ! number of BCs and dual grid node number
  integer              :: n_p, n_d, f_d, b_d, f_p, c_p, cnt
  integer, allocatable :: full_edge_n (:,:)   ! edges, nodes
  integer, allocatable :: full_edge_fb(:)     ! edges' faces at boundary
  integer, allocatable :: full_edge_bc(:)     ! edges' faces boundary colors
  integer, allocatable :: comp_edge_s(:)      ! compressed edge start
  integer, allocatable :: comp_edge_e(:)      ! compressed edge end
  integer, allocatable :: cell_to_node(:)
  integer, allocatable :: node_to_face(:)
  integer, allocatable :: node_to_cell(:)
  integer, allocatable :: edge_to_node(:)
  integer, allocatable :: edge_data(:)
  integer, allocatable :: cell_data(:)
  integer, allocatable :: node_data(:)
  integer, allocatable :: sharp_corner(:)     ! for sharp corners one day
  integer              :: c_p_list(2048)      ! prim cell and ...
  integer              :: n_d_list(2048)      ! ... dual node list
  integer              :: curr_f_d, curr_b_d, unused, dual_f_here
!==============================================================================!

  ! Alias(es)
  n_bc = prim % n_bnd_cond

  !-----------------------------!
  !                             !
  !   Set dual's problem_name   !
  !                             !
  !-----------------------------!
  problem_name(1) = trim(problem_name(1)) // '-dual'

  !-------------------------------------------------------------!
  !                                                             !
  !   Almost sure: you will be creating polyhedral grids here   !
  !                                                             !
  !-------------------------------------------------------------!
  dual % polyhedral = .true.

  !---------------------------------------!
  !                                       !
  !   Look for edges in the primal grid   !
  !                                       !
  !---------------------------------------!
  prim % n_edges = 0
  do s = 1, prim % n_faces
    prim % n_edges = prim % n_edges + prim % faces_n_nodes(s)
  end do
  print *, '# Expanded edges in primal grid:', prim % n_edges

  allocate(full_edge_n (prim % n_edges, 2));  full_edge_n (:,:) = 0
  allocate(full_edge_fb(prim % n_edges));     full_edge_fb(:)   = 0
  allocate(full_edge_bc(prim % n_edges));     full_edge_bc(:)   = 0

  !--------------------------------------------!
  !   Pile up all the edges on the prim grid   !
  !--------------------------------------------!
  prim % n_edges = 0
  do s = 1, prim % n_faces

    c1 = prim % faces_c(1, s)
    c2 = prim % faces_c(2, s)

    do i_nod = 1, prim % faces_n_nodes(s)
      j_nod = i_nod + 1;  if(j_nod > prim % faces_n_nodes(s)) j_nod = 1
      n1 = prim % faces_n(i_nod, s)
      n2 = prim % faces_n(j_nod, s)

      ! Increase the counter and store the data in it
      prim % n_edges = prim % n_edges + 1
      full_edge_n (prim % n_edges, 1) = min(n1, n2)
      full_edge_n (prim % n_edges, 2) = max(n1, n2)
      full_edge_fb(prim % n_edges)    =  s
      if(c2 > 0) then
        full_edge_bc(prim % n_edges) = 0  ! not a boundary face
      else
        full_edge_bc(prim % n_edges) = prim % bnd_cond % color(c2)  ! store color
      end if
    end do
  end do

  !--------------------------------------------------------!
  !   Sort all the edges by their node numbers and carry   !
  !   information on faces and boundary faces along        !
  !--------------------------------------------------------!
  call Sort_Mod_2_Int_Carry_2_Int(full_edge_n (1:prim % n_edges, 1),  &
                                  full_edge_n (1:prim % n_edges, 2),  &
                                  full_edge_fb(1:prim % n_edges),     &
                                  full_edge_bc(1:prim % n_edges))

  !-----------------------------------------!
  !   Estimate number of edges compressed   !
  !-----------------------------------------!
  prim % n_edges = 1
  do e = 2, size(full_edge_fb)  ! I don't like it, but here it is
    if(e > 1) then
      if(full_edge_n (e,1) .ne. full_edge_n (e-1,1) .or.  &
         full_edge_n (e,2) .ne. full_edge_n (e-1,2)) then
        prim % n_edges = prim % n_edges + 1
      end if
    end if
  end do
  print *, '# Compressed number of edges in primal grid:', prim % n_edges

  !----------------------------------------------------!
  !   Allocate memory for mapping and related arrays   !
  !----------------------------------------------------!
  allocate(comp_edge_s(prim % n_edges));  comp_edge_s(:) = 0
  allocate(comp_edge_e(prim % n_edges));  comp_edge_e(:) = 0

  allocate(prim % edges_n (2,      prim % n_edges))
  allocate(prim % edges_bc(0:n_bc, prim % n_edges))
  allocate(prim % edges_fb(2,      prim % n_edges))
  prim % edges_n (:,:) = 0
  prim % edges_bc(:,:) = 0  ! as if .false.
  prim % edges_fb(:,:) = 0

  !---------------------------------------------------------------!
  !   Store compressed edges with their information on boundary   !
  !---------------------------------------------------------------!
  prim % n_edges = 1
  comp_edge_s(prim % n_edges) = 1
  do e = 1, size(full_edge_fb)  ! I don't like it, but here it is

    if(e > 1) then
      if(full_edge_n (e, 1) .ne. full_edge_n (e-1, 1) .or.  &
         full_edge_n (e, 2) .ne. full_edge_n (e-1, 2)) then
        comp_edge_e(prim % n_edges) = e-1  ! mark the end of the old edge
        prim % n_edges = prim % n_edges + 1
        comp_edge_s(prim % n_edges) = e    ! mark the start of the new edge
      end if
    end if

    prim % edges_n (1, prim % n_edges) = full_edge_n(e, 1)
    prim % edges_n (2, prim % n_edges) = full_edge_n(e, 2)
    prim % edges_bc(full_edge_bc(e), prim % n_edges) = 1  ! as if .true.

    ! Saves boundary faces - needed to check if boundary
    ! face angles are too sharp at boundary edges
    if(full_edge_bc(e) .ne. 0) then
      if(prim % edges_fb(1,prim % n_edges) .eq. 0) then
         prim % edges_fb(1,prim % n_edges) = full_edge_fb(e)
      else if(prim % edges_fb(1,prim % n_edges) .ne. 0) then
        prim % edges_fb(2,prim % n_edges) = full_edge_fb(e)
      else
        print *, 'ERROR'
      end if
    end if

  end do
  comp_edge_e(prim % n_edges) = e-1  ! mark the end of the last edge

  !-----------------------!
  !                       !
  !   Allocate mappings   !
  !                       !
  !-----------------------!
  allocate(cell_to_node(-prim % n_bnd_cells  &
                        :prim % n_cells));  cell_to_node(:) = 0
  allocate(sharp_corner( prim % n_nodes));  sharp_corner(:) = 0
  allocate(node_to_face( prim % n_nodes));  node_to_face(:) = 0
  allocate(node_to_cell( prim % n_nodes));  node_to_cell(:) = 0
  allocate(edge_to_node( prim % n_edges));  edge_to_node(:) = 0
  allocate(cell_data(-prim % n_bnd_cells  &
                     :prim % n_cells));  cell_data(:) = 0
  allocate(edge_data( prim % n_edges));  edge_data(:) = 0
  allocate(node_data( prim % n_nodes));  node_data(:) = 0

  !----------------------!
  !   Plot sharp edges   !
  !----------------------!
  unused = N_Sharp_Edges(prim, edge_data)     ! find sharp edges
  call Save_Vtu_Edges(prim, edge_data)

  !-------------------------!
  !                         !
  !   Boundary conditions   !
  !    numbers and names    !
  !                         !
  !-------------------------!
  dual % n_bnd_cond = prim % n_bnd_cond
  allocate(dual % bnd_cond % name(dual % n_bnd_cond))
  do bc = 1, prim % n_bnd_cond
    dual % bnd_cond % name(bc) = prim % bnd_cond % name(bc)
    call To_Upper_Case(dual % bnd_cond % name(bc))
  end do

  !----------------------------------!
  !                                  !
  !   Inside mapping                 !
  !                                  !
  !   nodes(prim) =--> cells(dual)   !
  !   edges(prim) =--> faces(dual)   !
  !   cells(prim) =--> nodes(dual)   !
  !                                  !
  !----------------------------------!

  !----------------------------------------!
  !                                        !
  !   Allocate memory for inside mapping   !
  !                                        !
  !----------------------------------------!

  ! Count boundary cells (and boundary faces) for the dual grid
  ! (Remember that nodes in prim correspond to cells in dual)
  dual % n_bnd_cells = 0
  do bc = 1, prim % n_bnd_cond
    dual % n_bnd_cells = dual % n_bnd_cells  &
                       + N_Nodes_In_Bnd_Color(prim, bc, node_data)
  end do
  dual % n_faces = prim % n_edges  &   ! for faces inside
                 + dual % n_bnd_cells  ! for faces on the boundary
  dual % n_cells = prim % n_nodes
  dual % n_nodes = prim % n_cells                  &
                 + prim % n_bnd_cells              &
                 + N_Sharp_Edges(prim, edge_data)  &
                 + N_Sharp_Corners(prim, sharp_corner)

  call Allocate_Memory(dual)

  print *, '# Number of sharp corners = ', N_Sharp_Corners(prim, sharp_corner)

  !-----------------------------------------!
  !                                         !
  !   Browse through all the edges to map   !
  !                                         !
  !-----------------------------------------!
  d_nn = prim % n_cells      &
       + prim % n_bnd_cells

  do e = 1, prim % n_edges

    f_d = e

    ! Store cells surrounding each face in the dual grid ...
    ! ... which correspond to nodes in the prim grid
    n1 = prim % edges_n(1, e)
    n2 = prim % edges_n(2, e)
    dual % faces_c(1, f_d) = n1
    dual % faces_c(2, f_d) = n2

    ! Store nodes for each face in dual grid ...
    ! ... which are cells around each edge in prim
    cnt = 0
    do i_edg = comp_edge_s(e), comp_edge_e(e)
      f_p = full_edge_fb(i_edg)       ! get face index from the prim grid
      cnt = cnt + 1;  c_p_list(cnt) = prim % faces_c(1, f_p)
      cnt = cnt + 1;  c_p_list(cnt) = prim % faces_c(2, f_p)
    end do
    call Sort_Mod_Unique_Int(c_p_list(1:cnt), cnt)

    ! Transform prim cells to dual nodes.  There is a shift
    ! since cells go from negative values and nodes do not
    do i = 1, cnt
      if(c_p_list(i) < 0)  n_d_list(i) = c_p_list(i) + prim % n_bnd_cells + 1
      if(c_p_list(i) > 0)  n_d_list(i) = c_p_list(i) + prim % n_bnd_cells
    end do

    ! Form dual's faces' nodes
    dual % faces_n_nodes(f_d) = cnt
    dual % faces_n(1:cnt,f_d) = n_d_list(1:cnt)

    ! Copy node coordinates
    do i = 1, cnt
      c_p = c_p_list(i)  ! prim cell number
      n_d = n_d_list(i)  ! dual node number
      dual % xn(n_d) = prim % xc(c_p)
      dual % yn(n_d) = prim % yc(c_p)
      dual % zn(n_d) = prim % zc(c_p)
    end do

    ! Store cell to node mapping (prim cell to dual node)
    do i = 1, cnt
      c_p = c_p_list(i)  ! prim cell number
      n_d = n_d_list(i)  ! dual node number
      cell_to_node(c_p) = n_d
    end do

    !------------------------------------------------------------!
    !   If the face in a sharp corner, add one more node to it   !
    !------------------------------------------------------------!
    if(edge_data(e) .gt. 0) then

      ! Additional dual node number
      d_nn = d_nn + 1

      ! Add extra node to dual's faces' nodes
      dual % faces_n_nodes(f_d) = cnt + 1
      dual % faces_n(cnt + 1:cnt + 1, f_d) = d_nn

      ! Copy extra node coordinates
      n1 = prim % edges_n(1, e)
      n2 = prim % edges_n(2, e)
      dual % xn(d_nn) = (prim % xn(n1) + prim % xn(n2)) * 0.5  ! copy coords
      dual % yn(d_nn) = (prim % yn(n1) + prim % yn(n2)) * 0.5
      dual % zn(d_nn) = (prim % zn(n1) + prim % zn(n2)) * 0.5

      ! Store edge to node mapping (prim edge to dual node)
      edge_to_node(e) = d_nn

    end if

  end do  ! through edges

  !----------------------------------!
  !                                  !
  !   Boundary mapping               !
  !                                  !
  !----------------------------------!

  ! Update current number of dual faces -> equal to the number of faces inside
  curr_f_d = prim % n_edges
  curr_b_d = 0

  do bc = 1, prim % n_bnd_cond

    !----------------------------------------------------!
    !   Call this to mark boundary cells in this color   !
    !----------------------------------------------------!
    dual_f_here = N_Nodes_In_Bnd_Color(prim, bc, node_data)
    unused      = N_Bnd_Cells_In_Color(prim, bc, cell_data)
    unused      = N_Edges_On_Bnd_Color(prim, bc, edge_data)

    !-----------------------------------------!
    !   Find dual's boundary face, and dual   !
    !   boundary cell nodes from prim cells   !
    !-----------------------------------------!
    do c = -prim % n_bnd_cells, -1
      if(cell_data(c) .gt. 0) then

        ! Take the prim cell's nodes (these are from prim)
        do i_nod = 1, prim % cells_n_nodes(c)
          n_p = prim % cells_n(i_nod, c)

          ! Additional boundary face in the dual grid
          f_d  = curr_f_d + node_data(n_p)
          dual % faces_n_nodes(f_d) = dual % faces_n_nodes(f_d) + 1
          dual % faces_n(dual % faces_n_nodes(f_d), f_d) = cell_to_node(c)

          ! Additional boundary cell in the dual grid
          b_d  = curr_b_d - node_data(n_p)
          dual % cells_n_nodes(b_d) = dual % cells_n_nodes(b_d) + 1
          dual % cells_n(dual % cells_n_nodes(b_d), b_d) = cell_to_node(c)
          dual % bnd_cond % color(b_d) = bc

          ! Store node_to_face (for the next step, adding edges)
          node_to_face(n_p) = f_d
          node_to_cell(n_p) = b_d

          ! Store cells surrounding each face in the dual grid ...
          ! ... which correspond to nodes in the prim grid
          dual % faces_c(1, f_d) = n_p  ! link to cell inside
          dual % faces_c(2, f_d) = b_d  ! link to boundary cell
        end do

      end if
    end do

    ! Update current number of dual faces
    curr_f_d = curr_f_d + dual_f_here
    curr_b_d = curr_b_d - dual_f_here

    !---------------------------------------------------------!
    !   Add additional nodes to dual's face from prim edges   !
    !   Here we work on the faces already introduced above    !
    !---------------------------------------------------------!
    do e = 1, prim % n_edges
      if(edge_data(e) .gt. 0) then

        ! Take the prim edge's nodes (these are from prim)
        do i_nod = 1, 2
          n_p = prim % edges_n(i_nod, e)

          ! This node_to_face was stored in the previous step
          f_d = node_to_face(n_p)
          dual % faces_n_nodes(f_d) = dual % faces_n_nodes(f_d) + 1
          dual % faces_n(dual % faces_n_nodes(f_d), f_d) = edge_to_node(e)

          ! This node_to_cell was stored in the previous step
          b_d = node_to_cell(n_p)
          dual % cells_n_nodes(b_d) = dual % cells_n_nodes(b_d) + 1
          dual % cells_n(dual % cells_n_nodes(b_d), b_d) = edge_to_node(e)

          ! The grid has sharp corners, add them to boundary faces and cells
          if(sharp_corner(n_p) .gt. 0) then
            n_d = dual % n_nodes - sharp_corner(n_p) + 1
            dual % xn(n_d) = prim % xn(n_p)
            dual % yn(n_d) = prim % yn(n_p)
            dual % zn(n_d) = prim % zn(n_p)

            f_d = node_to_face(n_p)
            dual % faces_n_nodes(f_d) = dual % faces_n_nodes(f_d) + 1
            dual % faces_n(dual % faces_n_nodes(f_d), f_d) = n_d

            b_d = node_to_cell(n_p)
            dual % cells_n_nodes(b_d) = dual % cells_n_nodes(b_d) + 1
            dual % cells_n(dual % cells_n_nodes(b_d), b_d) = n_d
          end if
        end do
      end if
    end do

  end do

  !---------------------!
  !                     !
  !   Sort the points   !
  !                     !
  !---------------------!
  do s = 1, dual % n_faces
    call Sort_Face_Nodes(dual, s)
  end do

  !----------------------------------!
  !                                  !
  !   Store cells' faces and nodes   !
  !                                  !
  !----------------------------------!
  do s = 1, dual % n_faces
    c1 = dual % faces_c(1, s)
    c2 = dual % faces_c(2, s)

    !----------------------------------------------------!
    !   Store faces surrounding each cell in dual grid   !
    !----------------------------------------------------!
    dual % cells_n_faces(c1) = dual % cells_n_faces(c1) + 1
    dual % cells_n_faces(c2) = dual % cells_n_faces(c2) + 1
    dual % cells_f(dual % cells_n_faces(c1), c1) = s
    dual % cells_f(dual % cells_n_faces(c2), c2) = s

    !----------------------------------------------------!
    !   Store nodes surrounding each cell in dual grid   !
    !----------------------------------------------------!
    do i_nod = 1, dual % faces_n_nodes(s)
      n = dual % faces_n(i_nod, s)

      ! Handle cell 1
      do j_nod = 1, dual % cells_n_nodes(c1)
        n1 = dual % cells_n(j_nod, c1)
        if(n1 .eq. n) goto 1
      end do  ! j_nod
      dual % cells_n_nodes(c1) = dual % cells_n_nodes(c1) + 1
      dual % cells_n(dual % cells_n_nodes(c1), c1) = n
1     continue

      ! Handle cell 2
      do j_nod = 1, dual % cells_n_nodes(c2)
        n2 = dual % cells_n(j_nod, c2)
        if(n2 .eq. n) goto 2
      end do  ! j_nod
      dual % cells_n_nodes(c2) = dual % cells_n_nodes(c2) + 1
      dual % cells_n(dual % cells_n_nodes(c2), c2) = n
2     continue

    end do  ! do i_nod

  end do  ! do s

  !--------------------------------------------------------!
  !                                                        !
  !   Fix the number of nodes for polyhedral (all) cells   !
  !                                                        !
  !--------------------------------------------------------!
  do c = 1, dual % n_cells
    call Sort_Mod_Int(dual % cells_n(1:dual % cells_n_nodes(c),c))
    dual % cells_n_nodes(c) = -dual % cells_n_nodes(c)
  end do

  end subroutine
