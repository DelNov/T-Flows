!==============================================================================!
  subroutine Create_Dual(Convert, Prim, Dual)
!------------------------------------------------------------------------------!
!>  Creates a dual grid based on the given primal grid.  The dual grid concept
!>  is used here to create a polyhedral grid out of a tetrahedral.
!------------------------------------------------------------------------------!
!   Philosophy                                                                 !
!                                                                              !
!   * This routine converts a tetrahedral primal grid into a polyhedral dual   !
!     grid. The basic idea is that each primal node becomes a dual cell, each  !
!     primal edge becomes a dual face, and each primal cell becomes a dual     !
!     node.                                                                    !
!   * To make the dual grid suitable near boundaries and sharp geometric       !
!     features, additional dual nodes are introduced for boundary closure,     !
!     sharp edges, and sharp corners.                                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert     !! parent class
  type(Grid_Type)     :: Prim, Dual  !! primal and dual grid
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG = .false.
!-----------------------------------[Locals]-----------------------------------!
  integer              :: bc, e_p, c1_p, c2_p, c1_d, c2_d
  integer              :: n_p, n_d, n1_p, n2_p, n1_d, n2_d
  integer              :: i, i_nod, j_nod, i_edg, i_nor, nn
  integer              :: n_bc, d_nn  ! number of BCs and Dual grid node number
  integer              :: s_p, s_d, b_d, c_p, c_d, cnt
  integer              :: n_prim_sharp_c, n_prim_sharp_e
  integer, allocatable :: full_edge_n (:,:)  ! edges, nodes
  integer, allocatable :: full_edge_fb(:)    ! edges' faces at boundary
  integer, allocatable :: full_edge_bc(:)    ! edges' faces boundary regions
  integer, allocatable :: sorted_edge_f(:)   ! first and last occurrence in ...
  integer, allocatable :: sorted_edge_l(:)   ! ... in sorted full_edge_* list
  integer, allocatable :: cell_to_node(:)
  integer, allocatable :: node_to_face(:)
  integer, allocatable :: node_to_cell(:)
  integer, allocatable :: edge_to_node(:)
  integer, allocatable :: prim_bnd_cell_flag_in_reg(:)
  integer, allocatable :: prim_node_rank_in_reg(:)
  integer, allocatable :: prim_edge_flag_in_reg(:)
  integer, allocatable :: prim_sharp_node_rank(:)
  integer, allocatable :: prim_sharp_edge_flag(:)
  integer, allocatable :: sharp_inject(:)     ! sharp corner injected in Dual
  integer, allocatable :: concave_link(:,:)
  integer, allocatable :: concave_link_cnt(:)
  integer              :: c_p_list(2048)      ! prim cell and ...
  integer              :: n_d_list(2048)      ! ... dual node list
  integer              :: curr_s_d, curr_b_d, unused, dual_f_here
  logical              :: error_occured, found
  real                 :: nx, ny, nz, delta
  real,    allocatable :: prim_node_dx(:,:)
  real,    allocatable :: prim_node_dy(:,:)
  real,    allocatable :: prim_node_dz(:,:)
  real,    allocatable :: prim_node_d (:,:)
  integer, allocatable :: prim_node_n_norms(:)   ! number of normals
!==============================================================================!

  ! Alias(es)
  n_bc = Prim % n_bnd_regions

  !-----------------------------!
  !                             !
  !   Set Dual's problem_name   !
  !                             !
  !-----------------------------!
  problem_name(1) = trim(problem_name(1)) // '_dual'
  Dual % name = trim(Prim % name) // '_DUAL'

  !-------------------------------------------------------------!
  !                                                             !
  !   Almost sure: you will be creating polyhedral grids here   !
  !                                                             !
  !-------------------------------------------------------------!
  Dual % polyhedral = .true.

  !---------------------------------------!
  !                                       !
  !   Look for edges in the primal grid   !
  !                                       !
  !---------------------------------------!
  Prim % n_edges = 0
  do s_p = 1, Prim % n_faces
    Prim % n_edges = Prim % n_edges + Prim % faces_n_nodes(s_p)
  end do
  print *, '# Expanded edges in primal grid:', Prim % n_edges

  allocate(full_edge_n (Prim % n_edges, 2));  full_edge_n (:,:) = 0
  allocate(full_edge_fb(Prim % n_edges));     full_edge_fb(:)   = 0
  allocate(full_edge_bc(Prim % n_edges));     full_edge_bc(:)   = 0

  !--------------------------------------------!
  !   Pile up all the edges on the Prim grid   !
  !--------------------------------------------!
  Prim % n_edges = 0
  do s_p = 1, Prim % n_faces

    c1_p = Prim % faces_c(1, s_p)
    c2_p = Prim % faces_c(2, s_p)

    do i_nod = 1, Prim % faces_n_nodes(s_p)
      j_nod = i_nod + 1;  if(j_nod > Prim % faces_n_nodes(s_p)) j_nod = 1
      n1_p = Prim % faces_n(i_nod, s_p)
      n2_p = Prim % faces_n(j_nod, s_p)

      ! Increase the counter and store the data in it
      Prim % n_edges = Prim % n_edges + 1
      full_edge_n (Prim % n_edges, 1) = min(n1_p, n2_p)
      full_edge_n (Prim % n_edges, 2) = max(n1_p, n2_p)
      full_edge_fb(Prim % n_edges)    = s_p
      if(c2_p > 0) then
        ! Not a boundary face
        full_edge_bc(Prim % n_edges) = 0
      else
        ! Store region
        full_edge_bc(Prim % n_edges) = Prim % region % at_cell(c2_p)
      end if
    end do
  end do

  !--------------------------------------------------------!
  !   Sort all the edges by their node numbers and carry   !
  !   information on faces and boundary faces along        !
  !--------------------------------------------------------!
  call Sort % Two_Int_Carry_Two_Int(full_edge_n (1:Prim % n_edges, 1),  &
                                    full_edge_n (1:Prim % n_edges, 2),  &
                                    full_edge_fb(1:Prim % n_edges),     &
                                    full_edge_bc(1:Prim % n_edges))

  !-----------------------------------------!
  !   Estimate number of edges compressed   !
  !-----------------------------------------!
  Prim % n_edges = 1
  do e_p = 2, size(full_edge_fb)  ! I don't like it, but here it is
    if(e_p > 1) then
      if(full_edge_n(e_p, 1) .ne. full_edge_n(e_p - 1, 1) .or.  &
         full_edge_n(e_p, 2) .ne. full_edge_n(e_p - 1, 2)) then
        Prim % n_edges = Prim % n_edges + 1
      end if
    end if
  end do
  print *, '# Compressed number of edges in primal grid:', Prim % n_edges

  !----------------------------------------------------!
  !   Allocate memory for mapping and related arrays   !
  !----------------------------------------------------!
  allocate(sorted_edge_f(Prim % n_edges));  sorted_edge_f(:) = 0
  allocate(sorted_edge_l(Prim % n_edges));  sorted_edge_l(:) = 0

  allocate(Prim % edges_n (2,      Prim % n_edges))
  allocate(Prim % edges_bc(0:n_bc, Prim % n_edges))
  allocate(Prim % edges_fb(2,      Prim % n_edges))
  Prim % edges_n (:,:) = 0
  Prim % edges_bc(:,:) = 0  ! as if .false.
  Prim % edges_fb(:,:) = 0

  !---------------------------------------------------------------!
  !   Store compressed edges with their information on boundary   !
  !---------------------------------------------------------------!
  Prim % n_edges = 1
  sorted_edge_f(Prim % n_edges) = 1
  do e_p = 1, size(full_edge_fb)  ! I don't like it, but here it is

    if(e_p > 1) then
      if(full_edge_n (e_p, 1) .ne. full_edge_n (e_p-1, 1) .or.  &
         full_edge_n (e_p, 2) .ne. full_edge_n (e_p-1, 2)) then
        sorted_edge_l(Prim % n_edges) = e_p-1  ! mark the end of the old edge
        Prim % n_edges = Prim % n_edges + 1
        sorted_edge_f(Prim % n_edges) = e_p    ! mark the start of the new edge
      end if
    end if

    Prim % edges_n (1, Prim % n_edges) = full_edge_n(e_p, 1)
    Prim % edges_n (2, Prim % n_edges) = full_edge_n(e_p, 2)
    Prim % edges_bc(full_edge_bc(e_p), Prim % n_edges) = 1  ! as if .true.

    ! Saves boundary faces - needed to check if boundary
    ! face angles are too sharp at boundary edges
    if(full_edge_bc(e_p) .ne. 0) then
      if(Prim % edges_fb(1,Prim % n_edges) .eq. 0) then
        Prim % edges_fb(1, Prim % n_edges) = full_edge_fb(e_p)
      else if(Prim % edges_fb(1,Prim % n_edges) .ne. 0) then
        Prim % edges_fb(2, Prim % n_edges) = full_edge_fb(e_p)
      else
        print *, '# ERROR'
      end if
    end if

  end do
  sorted_edge_l(Prim % n_edges) = e_p-1  ! mark the end of the last edge

  ! De-allocate two full_edge_... arrays
  ! full_edge_n  -> not needed after compressed Prim % edges_n is built
  ! full_edge_bc -> not needed after Prim % edges_bc / Prim % edges_fb are built
  ! full_edge_fb -> still needed for mapping, for it contains list of faces
  !
  ! For example, after sorting, these arrays could have the values:
  !
  ! i_edg   full_edge_n(:,i_edg)   full_edge_bc(i_edg)   full_edge_fb(i_edg)
  ! ------------------------------------------------------------------------
  ! 21      edge (3,7)             4                     face 100
  ! 22      edge (3,7)             4                     face 124
  ! 23      edge (3,7)             4                     face 131
  ! 24      edge (3,7)             4                     face 155
  ! 25      edge (3,8)             5                     face 102
  deallocate(full_edge_n)
  deallocate(full_edge_bc)

  !-----------------------!
  !                       !
  !   Allocate mappings   !
  !                       !
  !-----------------------!
  allocate(cell_to_node(-Prim % n_bnd_cells  &
                        :Prim % n_cells));  cell_to_node(:) = 0
  allocate(node_to_face( Prim % n_nodes));  node_to_face(:) = 0
  allocate(node_to_cell( Prim % n_nodes));  node_to_cell(:) = 0
  allocate(edge_to_node( Prim % n_edges));  edge_to_node(:) = 0

  allocate(prim_sharp_edge_flag(Prim % n_edges))
  allocate(prim_sharp_node_rank( Prim % n_nodes))

  !--------------------------------!
  !   Find sharp edges and nodes   !
  !--------------------------------!
  n_prim_sharp_e = Convert % N_Sharp_Edges(Prim,  &
                                           prim_sharp_edge_flag)
  n_prim_sharp_c = Convert % N_Sharp_Corners(Prim,                  &
                                             prim_sharp_edge_flag,  &
                                             prim_sharp_node_rank)
  ! Plot sharp edges for checking
  call Prim % Save_Vtu_Edges(prim_sharp_edge_flag)  ! -1, 0, +1

  !-------------------------!
  !                         !
  !   Boundary conditions   !
  !    numbers and names    !
  !                         !
  !-------------------------!
  call Dual % Allocate_Regions(Prim % n_bnd_regions)
  do bc = 1, Prim % n_bnd_regions
    Dual % region % name(bc) = Prim % region % name(bc)
    call String % To_Upper_Case(Dual % region % name(bc))
  end do

  !----------------------------------!
  !                                  !
  !   Inside mapping                 !
  !                                  !
  !   nodes(Prim) =--> cells(Dual)   !
  !   edges(Prim) =--> faces(Dual)   !
  !   cells(Prim) =--> nodes(Dual)   !
  !                                  !
  !----------------------------------!

  !----------------------------------------!
  !                                        !
  !   Allocate memory for inside mapping   !
  !                                        !
  !----------------------------------------!

  ! Count boundary cells (and boundary faces) for the Dual grid
  ! (Remember that nodes in Prim correspond to cells in Dual)
  Dual % n_bnd_cells = 0

  ! Allocate memory to store node ranks in a region
  ! and the node ranks for sharp nodes (corners)
  allocate(prim_node_rank_in_reg(Prim % n_nodes))

  ! Increase total number of boundary cells in Dual
  do bc = 1, Prim % n_bnd_regions
    Dual % n_bnd_cells = Dual % n_bnd_cells                     &
                       + Convert % N_Nodes_In_Region(Prim, bc,  &
                                                     prim_node_rank_in_reg)
  end do

  ! Estimate number of faces, cells and nodes in the Dual grid ...
  Dual % n_faces = Prim % n_edges  &   ! for faces inside
                 + Dual % n_bnd_cells  ! for faces on the boundary
  Dual % n_cells = Prim % n_nodes
  Dual % n_nodes = Prim % n_cells       &
                 + Prim % n_bnd_cells   &
                 + n_prim_sharp_e       &
                 + n_prim_sharp_c

  ! ... and allocate memory for the Dual grid
  call Convert % Allocate_Memory(Dual)
  allocate(concave_link(6,  Dual % n_nodes));  concave_link(:,:)   = 0
  allocate(concave_link_cnt(Dual % n_nodes));  concave_link_cnt(:) = 0
  allocate(sharp_inject(    Dual % n_faces));  sharp_inject(:)     = 0

  ! Find sharp corners (nodes) in the Prim grid
  print *, '# Number of sharp corners = ', n_prim_sharp_c

  !----------------------------------------------!
  !                                              !
  !   Form internal Dual faces from Prim edges   !
  !                                              !
  !- - - - - - - - - - - - - - - - - - - - - - - +-------!
  !   * Each Prim edge becomes one internal Dual face.   !
  !   * The two end nodes of the Prim edge become the    !
  !     two neighbouring Dual cells of that face.        !
  !   * The Prim cells surrounding the edge become       !
  !     the nodes of the Dual face.                      !
  !------------------------------------------------------!
  d_nn = Prim % n_cells + Prim % n_bnd_cells  ! cells(Prim) =--> nodes(Dual)

  do e_p = 1, Prim % n_edges

    ! Face on the Dual grid (s_d) corresponds to the edge of the Prim grid
    s_d = e_p  ! edges(Prim) =--> faces(Dual)

    ! Store cells surrounding each face in the Dual grid ...
    ! ... which correspond to nodes in the Prim grid
    n1_p = Prim % edges_n(1, e_p)
    n2_p = Prim % edges_n(2, e_p)
    c1_d = n1_p  ! nodes(Prim) =--> cells(Dual)
    c2_d = n2_p  ! nodes(Prim) =--> cells(Dual)
    Dual % faces_c(1, s_d) = c1_d
    Dual % faces_c(2, s_d) = c2_d

    ! Store nodes for each face in Dual grid ...
    ! ... which are cells around each edge in Prim
    cnt = 0
    do i_edg = sorted_edge_f(e_p), sorted_edge_l(e_p)
      s_p = full_edge_fb(i_edg)       ! get face index from the Prim grid
      cnt = cnt + 1;  c_p_list(cnt) = Prim % faces_c(1, s_p)
      cnt = cnt + 1;  c_p_list(cnt) = Prim % faces_c(2, s_p)
    end do
    call Sort % Unique_Int(c_p_list(1:cnt), cnt)

    ! Transform Prim cells to Dual nodes.  There is a shift
    ! since cells go from negative values and nodes do not
    do i = 1, cnt
      if(c_p_list(i) < 0)  n_d_list(i) = c_p_list(i) + Prim % n_bnd_cells + 1
      if(c_p_list(i) > 0)  n_d_list(i) = c_p_list(i) + Prim % n_bnd_cells
    end do

    ! Form Dual's faces' nodes
    Dual % faces_n_nodes(s_d) = cnt
    call Enlarge % Matrix_Int(Dual % faces_n, i=(/1,cnt/))
    Dual % faces_n(1:cnt,s_d) = n_d_list(1:cnt)

    ! Copy node coordinates
    do i = 1, cnt
      c_p = c_p_list(i)  ! Prim cell number
      n_d = n_d_list(i)  ! Dual node number
      Dual % xn(n_d) = Prim % xc(c_p)
      Dual % yn(n_d) = Prim % yc(c_p)
      Dual % zn(n_d) = Prim % zc(c_p)
    end do

    ! Store cell to node mapping (Prim cell to Dual node)
    do i = 1, cnt
      c_p = c_p_list(i)  ! Prim cell number
      n_d = n_d_list(i)  ! Dual node number
      cell_to_node(c_p) = n_d
    end do

    !------------------------------------------------------------!
    !   If the face in a sharp corner, add one more node to it   !
    !------------------------------------------------------------!
    if(prim_sharp_edge_flag(e_p) .ne. 0) then  ! it is -1 or +1

      ! Additional Dual node number
      d_nn = d_nn + 1

      ! Add extra node to Dual's faces' nodes
      Dual % faces_n_nodes(s_d) = cnt + 1
      call Enlarge % Matrix_Int(Dual % faces_n, i=(/1,cnt+1/))
      Dual % faces_n(cnt + 1:cnt + 1, s_d) = d_nn

      ! Copy extra node coordinates
      n1_p = Prim % edges_n(1, e_p)
      n2_p = Prim % edges_n(2, e_p)

      ! Copy coordinates
      Dual % xn(d_nn) = (Prim % xn(n1_p) + Prim % xn(n2_p)) * 0.5
      Dual % yn(d_nn) = (Prim % yn(n1_p) + Prim % yn(n2_p)) * 0.5
      Dual % zn(d_nn) = (Prim % zn(n1_p) + Prim % zn(n2_p)) * 0.5

      ! For a concave sharp edge, the extra Dual node inserted at the
      ! edge midpoint can make the node ordering of the corresponding
      ! Dual face ambiguous.  Store the two boundary-side Dual nodes
      ! linked to this extra node; Sort_Face_Nodes uses this as a
      ! local sorting hint.
      !
      ! Note: tests showed that concave_link_cnt after this step is 1
      if(prim_sharp_edge_flag(e_p) .eq. -1) then
        ! print *, ' # Sharp edge:', e_p
        do i_edg = sorted_edge_f(e_p), sorted_edge_l(e_p)
          s_p = full_edge_fb(i_edg)       ! get face index from the Prim grid
          c2_p = Prim % faces_c(2, s_p)
          n2_d = c2_p + Prim % n_bnd_cells + 1
          if(c2_p .lt. 0) then
            cnt = concave_link_cnt(d_nn)

            ! If cnt==0, it checks pair (1,2), if cnt==1, it checks ...
            ! ... pair (3,4) and if cnt==2, it checks pair (5,6)
            if(     concave_link(cnt*2 + 1, d_nn) .eq. 0) then
                    concave_link(cnt*2 + 1, d_nn) = n2_d
            else if(concave_link(cnt*2 + 2, d_nn) .eq. 0) then
                    concave_link(cnt*2 + 2, d_nn) = n2_d
              cnt = cnt + 1
            end if

            concave_link_cnt(d_nn) = cnt
          end if
        end do
      end if

      ! Store edge to node mapping (Prim edge to Dual node)
      edge_to_node(e_p) = d_nn

    end if

  end do  ! through edges; e_p

  ! De-allocate arrays you won't need any more
  deallocate(sorted_edge_f)
  deallocate(sorted_edge_l)
  deallocate(full_edge_fb)

  !---------------------------------------!
  !   Sort the points on internal faces   !
  !---------------------------------------!
  do s_d = 1, Prim % n_edges  ! edges(Prim) =--> faces(Dual)
    call Convert % Sort_Face_Nodes(Dual, s_d, concave_link, concave_link_cnt)
  end do
  concave_link(:,:)   = 0  ! reset concave links since they ...
                           ! ... will now be used for inside faces
  concave_link_cnt(:) = 0

  !-------------------------------------------------------------------!
  !                                                                   !
  !   Form boundary Dual faces from boundary cells of the Prim grid   !
  !                                                                   !
  !-------------------------------------------------------------------!

  ! Update current number of Dual faces -> equal to the number of faces inside
  curr_s_d = Prim % n_edges
  curr_b_d = 0

  ! Allocate memory
  allocate(prim_edge_flag_in_reg(Prim % n_edges))
  prim_edge_flag_in_reg(:) = 0

  allocate(prim_bnd_cell_flag_in_reg(-Prim % n_bnd_cells:Prim % n_cells))
  prim_bnd_cell_flag_in_reg(:) = 0

  error_occured = .false.

  do bc = 1, Prim % n_bnd_regions

    !---------------------------------------------------------------!
    !   Call this to mark boundary nodes and cells in this region   !
    !---------------------------------------------------------------!
    dual_f_here = Convert % N_Nodes_In_Region    (Prim, bc,  &
                                                  prim_node_rank_in_reg)
    unused      = Convert % N_Bnd_Cells_In_Region(Prim, bc,  &
                                                  prim_bnd_cell_flag_in_reg)

    !-----------------------------------------!
    !   Find Dual's boundary face, and Dual   !
    !   boundary cell nodes from Prim cells   !
    !-----------------------------------------!
    do c_p = -Prim % n_bnd_cells, -1  ! boundary cells of the Prim grid
      if(prim_bnd_cell_flag_in_reg(c_p) .ne. 0) then

        ! Take the Prim cell's nodes (these are from Prim)
        do i_nod = 1, Prim % cells_n_nodes(c_p)
          n_p = Prim % cells_n(i_nod, c_p)

          ! Additional boundary face in the Dual grid
          s_d  = curr_s_d + prim_node_rank_in_reg(n_p)
          Dual % faces_n_nodes(s_d) = Dual % faces_n_nodes(s_d) + 1
          Dual % faces_n(Dual % faces_n_nodes(s_d), s_d) = cell_to_node(c_p)

          ! Additional boundary cell in the Dual grid
          b_d  = curr_b_d - prim_node_rank_in_reg(n_p)  ! boundary cells have
                                                        ! negative indices
          Dual % cells_n_nodes(b_d) = Dual % cells_n_nodes(b_d) + 1
          call Enlarge % Matrix_Int(Dual % cells_n,  &
                                    i=(/1,Dual % cells_n_nodes(b_d)/))
          Dual % cells_n(Dual % cells_n_nodes(b_d), b_d) = cell_to_node(c_p)
          Dual % region % at_cell(b_d) = bc

          ! Store node_to_face (for the next step, adding edges)
          node_to_face(n_p) = s_d
          node_to_cell(n_p) = b_d

          ! Store cells surrounding each face in the Dual grid ...
          ! ... which correspond to nodes in the Prim grid
          Dual % faces_c(1, s_d) = n_p  ! link to cell inside
          Dual % faces_c(2, s_d) = b_d  ! link to boundary cell
        end do

      end if
    end do

    ! Update current number of Dual faces and dual boundary cells
    curr_s_d = curr_s_d + dual_f_here
    curr_b_d = curr_b_d - dual_f_here

    !---------------------------------------------------------!
    !   Add additional nodes to Dual's face from Prim edges   !
    !      (Note that here we work on the faces already       !
    !       introduced above, just adding more nodes)         !
    !---------------------------------------------------------!
    unused = Convert % N_Edges_In_Region(Prim, bc, prim_edge_flag_in_reg)

    do e_p = 1, Prim % n_edges
      if(prim_edge_flag_in_reg(e_p) .ne. 0) then

        ! Take the Prim edge's nodes (these are from Prim)
        do i_nod = 1, 2
          n_p = Prim % edges_n(i_nod, e_p)

          ! This node_to_face was stored in the previous step
          s_d = node_to_face(n_p)
          Dual % faces_n_nodes(s_d) = Dual % faces_n_nodes(s_d) + 1
          call Enlarge % Matrix_Int(Dual % faces_n,  &
                                    i=(/1,Dual % faces_n_nodes(s_d)/))
          Dual % faces_n(Dual % faces_n_nodes(s_d), s_d) = edge_to_node(e_p)

          ! This node_to_cell was stored in the previous step
          b_d = node_to_cell(n_p)
          Dual % cells_n_nodes(b_d) = Dual % cells_n_nodes(b_d) + 1
          call Enlarge % Matrix_Int(Dual % cells_n,  &
                                    i=(/1,Dual % cells_n_nodes(b_d)/))
          Dual % cells_n(Dual % cells_n_nodes(b_d), b_d) = edge_to_node(e_p)
        end do  ! i_nod for edge, goes from 1 to 2

      end if
    end do  ! through edges

    !------------------------------------------------------------!
    !   Add additional nodes to Dual's face from sharp corners   !
    !        (Note that here we work on the faces already        !
    !         introduced above, just adding more nodes)          !
    !------------------------------------------------------------!
    do e_p = 1, Prim % n_edges
      if(prim_edge_flag_in_reg(e_p) .ne. 0) then

        do i_nod = 1, 2
          n_p = Prim % edges_n(i_nod, e_p)
          s_d = node_to_face(n_p)
          b_d = node_to_cell(n_p)

          Assert(s_d .gt. Prim % n_edges)

          ! The grid has sharp corners, add them to boundary faces and cells
          if(prim_sharp_node_rank(n_p) .gt. 0) then
            if(sharp_inject(s_d) .eq. 0) then

              ! Estimate the number of new node in Dual
              ! (prim_sharp_node_rank holds its local number)
              n_d = Dual % n_nodes - prim_sharp_node_rank(n_p) + 1
              Dual % xn(n_d) = Prim % xn(n_p)
              Dual % yn(n_d) = Prim % yn(n_p)
              Dual % zn(n_d) = Prim % zn(n_p)

              ! Insert a new node from the edge to the bundary face
              Dual % faces_n_nodes(s_d) = Dual % faces_n_nodes(s_d) + 1
              call Enlarge % Matrix_Int(Dual % faces_n,  &
                                        i=(/1,Dual % faces_n_nodes(s_d)/))
              Dual % faces_n(Dual % faces_n_nodes(s_d), s_d) = n_d

              ! Insert a new node from the edge to the bundary cell
              Dual % cells_n_nodes(b_d) = Dual % cells_n_nodes(b_d) + 1
              call Enlarge % Matrix_Int(Dual % cells_n,  &
                                        i=(/1,Dual % cells_n_nodes(b_d)/))
              Dual % cells_n(Dual % cells_n_nodes(b_d), b_d) = n_d

              ! Mark that the face has been injected a sharp corner
              ! write(*,'(a,i9)') ' # Injecting new node in face', s_d
              sharp_inject(s_d) = sharp_inject(s_d) + 1

              ! Check those little nodes inserted just before
              n1_d = Dual % faces_n(Dual % faces_n_nodes(s_d)-2, s_d)
              n2_d = Dual % faces_n(Dual % faces_n_nodes(s_d)-1, s_d)

              ! Increase the concave link count for this node
              concave_link_cnt(n_d) = concave_link_cnt(n_d) + 1

              ! Store the first concave link for this node
              cnt = concave_link_cnt(n_d)

              ! Fill paris 1,2 or 3,4 or 5,6
              concave_link(cnt*2 - 1, n_d) = n1_d
              concave_link(cnt*2,     n_d) = n2_d

              ! A sharp-corner Dual node can be reached from more than one
              ! boundary edge. This is expected at domain vertices where
              ! several sharp edges meet, for example at the eight corners
              ! of a cube.  In that case concave_link_cnt(n_d) counts how many
              ! edge-neighbour pairs were attached to the same sharp-corner
              ! node.  The pairs are stored in concave_link at first index
              ! positions (1,2), (3,4), and (5,6).
              !
              ! Note: tests showed that concave_link_cnt after this step is 3
              if(concave_link_cnt(n_d) .gt. 3) then
                error_occured = .true.
              end if
            end if  ! sharp inject in s_d
          end if    ! sharp corner here

        end do  ! i_nod for edge, goes from 1 to 2
      end if    ! edge is in the region
    end do      ! through edges

  end do  ! bc, boundary region

  if(error_occured) then
    call Message % Error(66, "Dual nodes which are sharp corners "     //     &
                             "(such a eight cube corners) with more "  //     &
                             "than three coinciding edges found. "     //     &
                             "\n \n "                                  //     &
                             "Such a situation haven't been predicted "//     &
                             "and the program will have to stop. "     //     &
                             "Contact the authors and send them this grid.",  &
                             file = __FILE__, line = __LINE__)
  end if

  ! De-allocate what you won't need any more
  deallocate(sharp_inject)
  deallocate(prim_sharp_node_rank)
  deallocate(prim_node_rank_in_reg)
  deallocate(edge_to_node)
  deallocate(node_to_cell)
  deallocate(node_to_face)
  deallocate(cell_to_node)

  !---------------------------------------!
  !                                       !
  !   Sort the points on boundary faces   !
  !                                       !
  !---------------------------------------!
  do s_d = Prim % n_edges + 1, Dual % n_faces
    call Convert % Sort_Face_Nodes(Dual, s_d, concave_link, concave_link_cnt)
  end do

  deallocate(concave_link_cnt)
  deallocate(concave_link)

  !------------------------------------------------------------!
  !   Save only internal faces to see if sorting of boundary   !
  !       faces messed up already sorted internal faces        !
  !------------------------------------------------------------!
  ! Dual % n_faces = Prim % n_edges
  ! call Dual % Save_Vtu_Faces()
  ! stop

  !---------------------------------------------------!
  !   Save all faces to see if sorting of boundary    !
  !   faces messed up already sorted internal faces   !
  !---------------------------------------------------!
  ! call Dual % Save_Vtu_Faces()
  ! stop

  !----------------------------------!
  !                                  !
  !   Store cells' faces and nodes   !
  !                                  !
  !----------------------------------!
  do s_d = 1, Dual % n_faces
    c1_d = Dual % faces_c(1, s_d)
    c2_d = Dual % faces_c(2, s_d)

    !----------------------------------------------------!
    !   Store faces surrounding each cell in Dual grid   !
    !----------------------------------------------------!
    Dual % cells_n_faces(c1_d) = Dual % cells_n_faces(c1_d) + 1
    Dual % cells_n_faces(c2_d) = Dual % cells_n_faces(c2_d) + 1
    call Enlarge % Matrix_int(Dual % cells_f,i=(/1,Dual % cells_n_faces(c1_d)/))
    call Enlarge % Matrix_int(Dual % cells_f,i=(/1,Dual % cells_n_faces(c2_d)/))
    Dual % cells_f(Dual % cells_n_faces(c1_d), c1_d) = s_d
    Dual % cells_f(Dual % cells_n_faces(c2_d), c2_d) = s_d

    !----------------------------------------------------!
    !   Store nodes surrounding each cell in Dual grid   !
    !----------------------------------------------------!
    do i_nod = 1, Dual % faces_n_nodes(s_d)
      n_d = Dual % faces_n(i_nod, s_d)

      ! Handle cell 1
      do j_nod = 1, Dual % cells_n_nodes(c1_d)
        n1_d = Dual % cells_n(j_nod, c1_d)
        if(n1_d .eq. n_d) goto 1
      end do  ! j_nod
      Dual % cells_n_nodes(c1_d) = Dual % cells_n_nodes(c1_d) + 1
      call Enlarge % Matrix_Int(Dual % cells_n,  &
                                i=(/1,Dual % cells_n_nodes(c1_d)/))
      Dual % cells_n(Dual % cells_n_nodes(c1_d), c1_d) = n_d
1     continue

      ! Handle cell 2
      do j_nod = 1, Dual % cells_n_nodes(c2_d)
        n2_d = Dual % cells_n(j_nod, c2_d)
        if(n2_d .eq. n_d) goto 2
      end do  ! j_nod
      Dual % cells_n_nodes(c2_d) = Dual % cells_n_nodes(c2_d) + 1
      call Enlarge % Matrix_Int(Dual % cells_n,  &
                                i=(/1,Dual % cells_n_nodes(c2_d)/))
      Dual % cells_n(Dual % cells_n_nodes(c2_d), c2_d) = n_d
2     continue

    end do  ! do i_nod

  end do  ! do s

  !--------------------------------------------------------!
  !                                                        !
  !   Fix the number of nodes for polyhedral (all) cells   !
  !                                                        !
  !--------------------------------------------------------!
  do c_d = 1, Dual % n_cells
    call Sort % Int_Array(Dual % cells_n(1:Dual % cells_n_nodes(c_d),c_d))
    Dual % cells_n_nodes(c_d) = -Dual % cells_n_nodes(c_d)
    Assert(Dual % cells_n_nodes(c_d) .lt. 0)
  end do

  !------------------------!
  !                        !
  !   Mark concave cells   !
  !                        !
  !- - - - - - - - - - - - +---------------------------------------------!
  !   Please note that Prim node -> Dual cell works directly because I   !
  !   preserved the numbering.                                           !
  !                                                                      !
  !   However, Prim cell ->  Dual node would not work directly because   !
  !   dual nodes are shifted by boundary nodes and mixed with extra      !
  !   sharp-feature nodes.  (Think about it a bit more and try to        !
  !   explain it in simpler words.)                                      !
  !----------------------------------------------------------------------!
  do e_p = 1, Prim % n_edges
    if(prim_sharp_edge_flag(e_p) .lt. 0) then
      c1_d = Prim % edges_n(1, e_p)   ! nodes(Prim) =--> cells(Dual)
      c2_d = Prim % edges_n(2, e_p)   ! nodes(Prim) =--> cells(Dual)
      Assert(c1_d .le. Dual % n_cells)
      Assert(c2_d .le. Dual % n_cells)
      Dual % concave(c1_d) = .true.
      Dual % concave(c2_d) = .true.
    end if
  end do

  ! Allocate memory for helping array
  allocate(prim_node_n_norms(Prim % n_nodes));  prim_node_n_norms(:) = 0

  !-----------------------------------------------------------------------!
  !                                                                       !
  !   For each "concave" Prim node, find characterstic boundary normals   !
  !                                                                       !
  !-----------------------------------------------------------------------!
  do s_p = 1, Prim % n_faces
    c1_p = Prim % faces_c(1, s_p)
    c2_p = Prim % faces_c(2, s_p)
    if(c2_p .lt. 0) then
      do i_nod = 1, Prim % faces_n_nodes(s_p)
        n_p = Prim % faces_n(i_nod, s_p)
        c_d = n_p  ! nodes(Prim) =--> cells(Dual)
        if(Dual % concave(c_d)) then

          ! Surface magnitude is not yet computed at this point
          Prim % s(s_p) = sqrt(  Prim % sx(s_p) ** 2  &
                               + Prim % sy(s_p) ** 2  &
                               + Prim % sz(s_p) ** 2)
          nx = Prim % sx(s_p) / Prim % s(s_p)
          ny = Prim % sy(s_p) / Prim % s(s_p)
          nz = Prim % sz(s_p) / Prim % s(s_p)

          nn = prim_node_n_norms(n_p)

          !---------------------------------------------!
          !   If no normals have been found/saved yet   !
          !---------------------------------------------!
          if(nn .eq. 0) then

            prim_node_n_norms(n_p) = 1  ! first normal
            call Enlarge % Matrix_Real(prim_node_dx,  &
                                       i=(/1,1/),     &
                                       j=(/1,Prim % n_nodes/))
            call Enlarge % Matrix_Real(prim_node_dy,  &
                                       i=(/1,1/),     &
                                       j=(/1,Prim % n_nodes/))
            call Enlarge % Matrix_Real(prim_node_dz,  &
                                       i=(/1,1/),     &
                                       j=(/1,Prim % n_nodes/))
            call Enlarge % Matrix_Real(prim_node_d,   &
                                       i=(/1,1/),     &
                                       j=(/1,Prim % n_nodes/))
            prim_node_dx(1, n_p) = nx
            prim_node_dy(1, n_p) = ny
            prim_node_dz(1, n_p) = nz
            prim_node_d (1, n_p) = sqrt(Prim % s(s_p))  ! measure of length

          !----------------------------------------------------!
          !   Some normals have been already found and saved   !
          !----------------------------------------------------!
          else

            ! Check if this normal has been stored already
            found = .false.
            do i_nor = 1, nn  ! browse through stored normals
              if(Math % Approx_Real(nx, prim_node_dx(i_nor, n_p), MICRO) .and. &
                 Math % Approx_Real(ny, prim_node_dy(i_nor, n_p), MICRO) .and. &
                 Math % Approx_Real(nz, prim_node_dz(i_nor, n_p), MICRO)) then
                found = .true.
              end if
            end do

            ! If not found, store one more normal
            if(.not. found) then
              prim_node_n_norms(n_p) = prim_node_n_norms(n_p) + 1
              nn = prim_node_n_norms(n_p)
              call Enlarge % Matrix_Real(prim_node_dx,  &
                                         i=(/1,nn/),    &
                                         j=(/1,Prim % n_nodes/))
              call Enlarge % Matrix_Real(prim_node_dy,  &
                                         i=(/1,nn/),    &
                                         j=(/1,Prim % n_nodes/))
              call Enlarge % Matrix_Real(prim_node_dz,  &
                                         i=(/1,nn/),    &
                                         j=(/1,Prim % n_nodes/))
              call Enlarge % Matrix_Real(prim_node_d,   &
                                         i=(/1,nn/),    &
                                         j=(/1,Prim % n_nodes/))
              prim_node_dx(nn, n_p) = nx
              prim_node_dy(nn, n_p) = ny
              prim_node_dz(nn, n_p) = nz
              prim_node_d (nn, n_p) = sqrt(Prim % s(s_p))  ! measure of length
            end if
          end if
        end if
      end do
    end if
  end do

  do c_d = 1, Dual % n_cells
    if(Dual % concave(c_d)) then
      n_p = c_d  ! nodes(Prim) =--> cells(Dual)

      ! First find characteristic shift
      delta = 0
      do i_nor = 1, prim_node_n_norms(n_p)
        delta = delta + prim_node_d(i_nor, n_p)
      end do
      delta = delta / prim_node_n_norms(n_p)
      delta = delta * 0.2  ! the 0.2 is an ad-hoc guess

      ! Just copy the coordinates of the Prim node ...
      Dual % xc(c_d) = Prim % xn(n_p)
      Dual % yc(c_d) = Prim % yn(n_p)
      Dual % zc(c_d) = Prim % zn(n_p)

      ! .. then shift it according to normals
      do i_nor = 1, prim_node_n_norms(n_p)
        Dual % xc(c_d) = Dual % xc(c_d) - prim_node_dx(i_nor, n_p) * delta
        Dual % yc(c_d) = Dual % yc(c_d) - prim_node_dy(i_nor, n_p) * delta
        Dual % zc(c_d) = Dual % zc(c_d) - prim_node_dz(i_nor, n_p) * delta
      end do
      if(DEBUG) then
        call Dual % Save_Vtk_Cell(c_d, "concave_cell", c_d, plot_center=.true.)
      end if
    end if
  end do

  end subroutine
