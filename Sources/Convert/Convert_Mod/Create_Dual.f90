!==============================================================================!
  subroutine Create_Dual(Convert, Prim, Dual)
!------------------------------------------------------------------------------!
!>  Creates a dual grid based on the given primal grid.  The dual grid concept
!>  is used here to create a polyhedral grid out of a tetrahedral.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initial setup: Sets up the problem name for the dual grid and marks it   !
!     as polyhedral.                                                           !
!   * Edge processing: Looks for edges in the primal grid and stores their     !
!     information.                                                             !
!   * Edge sorting and compression: Sorts all primal edges by their node       !
!     numbers and compresses them to eliminate duplicates.                     !
!   * Allocate mappings and arrays: Allocates memory for various arrays and    !
!     mappings, such as those relating edges to nodes and cells to nodes.      !
!   * Boundary conditions processing: Processes boundary conditions, including !
!     updating the number of boundary cells and faces in the dual grid.        !
!   * Inside mapping: Maps nodes of the primal grid to cells of the dual grid, !
!     edges of the primal to faces of the dual, and cells of the primal to     !
!     nodes of the dual.                                                       !
!   * Memory allocation for dual grid: Allocates memory for the dual grid      !
!     based on the calculated sizes of various components.                     !
!   * Handling sharp edges and corners: Processes sharp edges and corners in   !
!     the primal grid and reflects these in the dual grid structure.           !
!   * Face and cell assembly in dual grid: Assembles faces and cells in the    !
!     dual grid based on the primal-to-dual mappings.                          !
!   * Final Adjustments: Makes final adjustments to the structure of the dual  !
!     grid, ensuring that it accurately represents the topology of the primal  !
!     grid.                                                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert     !! parent class
  type(Grid_Type)     :: Prim, Dual  !! primal and dual grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: bc, e, c, c1, c2, s, n, n1, n2
  integer              :: i, i_nod, j_nod, i_edg
  integer              :: n_bc, d_nn  ! number of BCs and Dual grid node number
  integer              :: n_p, n_d, f_d, b_d, f_p, c_p, cnt
  integer, allocatable :: full_edge_n (:,:)   ! edges, nodes
  integer, allocatable :: full_edge_fb(:)     ! edges' faces at boundary
  integer, allocatable :: full_edge_bc(:)     ! edges' faces boundary regions
  integer, allocatable :: comp_edge_f(:)      ! compressed edge first
  integer, allocatable :: comp_edge_l(:)      ! compressed edge last
  integer, allocatable :: cell_to_node(:)
  integer, allocatable :: node_to_face(:)
  integer, allocatable :: node_to_cell(:)
  integer, allocatable :: edge_to_node(:)
  integer, allocatable :: edge_data(:)
  integer, allocatable :: cell_data(:)
  integer, allocatable :: node_data(:)
  integer, allocatable :: sharp_corner(:)     ! sharp corners in Prim
  integer, allocatable :: sharp_inject(:)     ! sharp corner injected in Dual
  integer, allocatable :: concave_link(:,:)
  integer              :: c_p_list(2048)      ! Prim cell and ...
  integer              :: n_d_list(2048)      ! ... Dual node list
  integer              :: curr_f_d, curr_b_d, unused, dual_f_here
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
  do s = 1, Prim % n_faces
    Prim % n_edges = Prim % n_edges + Prim % faces_n_nodes(s)
  end do
  print *, '# Expanded edges in primal grid:', Prim % n_edges

  allocate(full_edge_n (Prim % n_edges, 2));  full_edge_n (:,:) = 0
  allocate(full_edge_fb(Prim % n_edges));     full_edge_fb(:)   = 0
  allocate(full_edge_bc(Prim % n_edges));     full_edge_bc(:)   = 0

  !--------------------------------------------!
  !   Pile up all the edges on the Prim grid   !
  !--------------------------------------------!
  Prim % n_edges = 0
  do s = 1, Prim % n_faces

    c1 = Prim % faces_c(1, s)
    c2 = Prim % faces_c(2, s)

    do i_nod = 1, Prim % faces_n_nodes(s)
      j_nod = i_nod + 1;  if(j_nod > Prim % faces_n_nodes(s)) j_nod = 1
      n1 = Prim % faces_n(i_nod, s)
      n2 = Prim % faces_n(j_nod, s)

      ! Increase the counter and store the data in it
      Prim % n_edges = Prim % n_edges + 1
      full_edge_n (Prim % n_edges, 1) = min(n1, n2)
      full_edge_n (Prim % n_edges, 2) = max(n1, n2)
      full_edge_fb(Prim % n_edges)    =  s
      if(c2 > 0) then
        full_edge_bc(Prim % n_edges) = 0  ! not a boundary face
      else
        full_edge_bc(Prim % n_edges) = Prim % region % at_cell(c2)  ! store reg
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
  do e = 2, size(full_edge_fb)  ! I don't like it, but here it is
    if(e > 1) then
      if(full_edge_n (e,1) .ne. full_edge_n (e-1,1) .or.  &
         full_edge_n (e,2) .ne. full_edge_n (e-1,2)) then
        Prim % n_edges = Prim % n_edges + 1
      end if
    end if
  end do
  print *, '# Compressed number of edges in primal grid:', Prim % n_edges

  !----------------------------------------------------!
  !   Allocate memory for mapping and related arrays   !
  !----------------------------------------------------!
  allocate(comp_edge_f(Prim % n_edges));  comp_edge_f(:) = 0
  allocate(comp_edge_l(Prim % n_edges));  comp_edge_l(:) = 0

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
  comp_edge_f(Prim % n_edges) = 1
  do e = 1, size(full_edge_fb)  ! I don't like it, but here it is

    if(e > 1) then
      if(full_edge_n (e, 1) .ne. full_edge_n (e-1, 1) .or.  &
         full_edge_n (e, 2) .ne. full_edge_n (e-1, 2)) then
        comp_edge_l(Prim % n_edges) = e-1  ! mark the end of the old edge
        Prim % n_edges = Prim % n_edges + 1
        comp_edge_f(Prim % n_edges) = e    ! mark the start of the new edge
      end if
    end if

    Prim % edges_n (1, Prim % n_edges) = full_edge_n(e, 1)
    Prim % edges_n (2, Prim % n_edges) = full_edge_n(e, 2)
    Prim % edges_bc(full_edge_bc(e), Prim % n_edges) = 1  ! as if .true.

    ! Saves boundary faces - needed to check if boundary
    ! face angles are too sharp at boundary edges
    if(full_edge_bc(e) .ne. 0) then
      if(Prim % edges_fb(1,Prim % n_edges) .eq. 0) then
         Prim % edges_fb(1,Prim % n_edges) = full_edge_fb(e)
      else if(Prim % edges_fb(1,Prim % n_edges) .ne. 0) then
        Prim % edges_fb(2,Prim % n_edges) = full_edge_fb(e)
      else
        print *, '# ERROR'
      end if
    end if

  end do
  comp_edge_l(Prim % n_edges) = e-1  ! mark the end of the last edge

  !-----------------------!
  !                       !
  !   Allocate mappings   !
  !                       !
  !-----------------------!
  allocate(cell_to_node(-Prim % n_bnd_cells  &
                        :Prim % n_cells));  cell_to_node(:) = 0
  allocate(sharp_corner( Prim % n_nodes));  sharp_corner(:) = 0
  allocate(node_to_face( Prim % n_nodes));  node_to_face(:) = 0
  allocate(node_to_cell( Prim % n_nodes));  node_to_cell(:) = 0
  allocate(edge_to_node( Prim % n_edges));  edge_to_node(:) = 0
  allocate(cell_data(-Prim % n_bnd_cells  &
                     :Prim % n_cells));  cell_data(:) = 0
  allocate(edge_data( Prim % n_edges));  edge_data(:) = 0
  allocate(node_data( Prim % n_nodes));  node_data(:) = 0

  !----------------------!
  !   Plot sharp edges   !
  !----------------------!
  unused = Convert % N_Sharp_Edges(Prim, edge_data)     ! find sharp edges
  call Prim % Save_Vtu_Edges(edge_data)

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
  do bc = 1, Prim % n_bnd_regions
    Dual % n_bnd_cells = Dual % n_bnd_cells  &
                       + Convert % N_Nodes_In_Region(Prim, bc, node_data)
  end do
  Dual % n_faces = Prim % n_edges  &   ! for faces inside
                 + Dual % n_bnd_cells  ! for faces on the boundary
  Dual % n_cells = Prim % n_nodes
  Dual % n_nodes = Prim % n_cells                  &
                 + Prim % n_bnd_cells              &
                 + Convert % N_Sharp_Edges(Prim, edge_data)  &
                 + Convert % N_Sharp_Corners(Prim, sharp_corner)

  call Convert % Allocate_Memory(Dual)
  allocate(concave_link(2, Dual % n_nodes));  concave_link(:,:) = 0
  allocate(sharp_inject(   Dual % n_faces));  sharp_inject(:)   = 0

  print *, '# Number of sharp corners = ',  &
            Convert % N_Sharp_Corners(Prim, sharp_corner)

  !-----------------------------------------!
  !                                         !
  !   Browse through all the edges to map   !
  !                                         !
  !-----------------------------------------!
  d_nn = Prim % n_cells      &
       + Prim % n_bnd_cells

  do e = 1, Prim % n_edges

    f_d = e

    ! Store cells surrounding each face in the Dual grid ...
    ! ... which correspond to nodes in the Prim grid
    n1 = Prim % edges_n(1, e)
    n2 = Prim % edges_n(2, e)
    Dual % faces_c(1, f_d) = n1
    Dual % faces_c(2, f_d) = n2

    ! Store nodes for each face in Dual grid ...
    ! ... which are cells around each edge in Prim
    cnt = 0
    do i_edg = comp_edge_f(e), comp_edge_l(e)
      f_p = full_edge_fb(i_edg)       ! get face index from the Prim grid
      cnt = cnt + 1;  c_p_list(cnt) = Prim % faces_c(1, f_p)
      cnt = cnt + 1;  c_p_list(cnt) = Prim % faces_c(2, f_p)
    end do
    call Sort % Unique_Int(c_p_list(1:cnt), cnt)

    ! Transform Prim cells to Dual nodes.  There is a shift
    ! since cells go from negative values and nodes do not
    do i = 1, cnt
      if(c_p_list(i) < 0)  n_d_list(i) = c_p_list(i) + Prim % n_bnd_cells + 1
      if(c_p_list(i) > 0)  n_d_list(i) = c_p_list(i) + Prim % n_bnd_cells
    end do

    ! Form Dual's faces' nodes
    Dual % faces_n_nodes(f_d) = cnt
    call Adjust_First_Dim(cnt, Dual % faces_n)
    Dual % faces_n(1:cnt,f_d) = n_d_list(1:cnt)

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
    if(edge_data(e) .ne. 0) then

      ! Additional Dual node number
      d_nn = d_nn + 1

      ! Add extra node to Dual's faces' nodes
      Dual % faces_n_nodes(f_d) = cnt + 1
      call Adjust_First_Dim(cnt + 1, Dual % faces_n)
      Dual % faces_n(cnt + 1:cnt + 1, f_d) = d_nn

      ! Copy extra node coordinates
      n1 = Prim % edges_n(1, e)
      n2 = Prim % edges_n(2, e)
      Dual % xn(d_nn) = (Prim % xn(n1) + Prim % xn(n2)) * 0.5  ! copy coords
      Dual % yn(d_nn) = (Prim % yn(n1) + Prim % yn(n2)) * 0.5
      Dual % zn(d_nn) = (Prim % zn(n1) + Prim % zn(n2)) * 0.5

      ! Mark node in the concave corner.  It seems impossible
      ! sort out the direction of nodes for a concave face
      ! Retreive boundary face information again
      if(edge_data(e) .eq. -1) then
        ! print *, ' # Sharp edge:', e
        do i_edg = comp_edge_f(e), comp_edge_l(e)
          f_p = full_edge_fb(i_edg)       ! get face index from the Prim grid
          c2 = Prim % faces_c(2, f_p)
          n2 = c2 + Prim % n_bnd_cells + 1
          if(c2 < 0) then
            if(concave_link(1,d_nn) .eq. 0) then
              concave_link(1,d_nn) = n2
            else if(concave_link(2,d_nn) .eq. 0) then
              concave_link(2,d_nn) = n2
            else
              print *, '# Well, well, this is pretty wrong!'
            end if
          end if
        end do
      end if

      ! Store edge to node mapping (Prim edge to Dual node)
      edge_to_node(e) = d_nn

    end if

  end do  ! through edges

  !---------------------------------------!
  !                                       !
  !   Sort the points on internal faces   !
  !                                       !
  !---------------------------------------!
  do s = 1, Prim % n_edges
    call Convert % Sort_Face_Nodes(Dual, s, concave_link)
  end do
  concave_link(:,:) = 0  ! reset concave links since they work only for inside faces

  !----------------------------------!
  !                                  !
  !   Boundary mapping               !
  !                                  !
  !----------------------------------!

  ! Update current number of Dual faces -> equal to the number of faces inside
  curr_f_d = Prim % n_edges
  curr_b_d = 0

  do bc = 1, Prim % n_bnd_regions

    !-----------------------------------------------------!
    !   Call this to mark boundary cells in this region   !
    !-----------------------------------------------------!
    dual_f_here = Convert % N_Nodes_In_Region(Prim, bc, node_data)
    unused      = Convert % N_Bnd_Cells_In_Region(Prim, bc, cell_data)
    unused      = Convert % N_Edges_In_Region(Prim, bc, edge_data)

    !-----------------------------------------!
    !   Find Dual's boundary face, and Dual   !
    !   boundary cell nodes from Prim cells   !
    !-----------------------------------------!
    do c = -Prim % n_bnd_cells, -1
      if(cell_data(c) .ne. 0) then

        ! Take the Prim cell's nodes (these are from Prim)
        do i_nod = 1, Prim % cells_n_nodes(c)
          n_p = Prim % cells_n(i_nod, c)

          ! Additional boundary face in the Dual grid
          f_d  = curr_f_d + node_data(n_p)
          Dual % faces_n_nodes(f_d) = Dual % faces_n_nodes(f_d) + 1
          Dual % faces_n(Dual % faces_n_nodes(f_d), f_d) = cell_to_node(c)

          ! Additional boundary cell in the Dual grid
          b_d  = curr_b_d - node_data(n_p)
          Dual % cells_n_nodes(b_d) = Dual % cells_n_nodes(b_d) + 1
          call Adjust_First_Dim(Dual % cells_n_nodes(b_d), Dual % cells_n)
          Dual % cells_n(Dual % cells_n_nodes(b_d), b_d) = cell_to_node(c)
          Dual % region % at_cell(b_d) = bc

          ! Store node_to_face (for the next step, adding edges)
          node_to_face(n_p) = f_d
          node_to_cell(n_p) = b_d

          ! Store cells surrounding each face in the Dual grid ...
          ! ... which correspond to nodes in the Prim grid
          Dual % faces_c(1, f_d) = n_p  ! link to cell inside
          Dual % faces_c(2, f_d) = b_d  ! link to boundary cell
        end do

      end if
    end do

    ! Update current number of Dual faces
    curr_f_d = curr_f_d + dual_f_here
    curr_b_d = curr_b_d - dual_f_here

    !---------------------------------------------------------!
    !   Add additional nodes to Dual's face from Prim edges   !
    !   Here we work on the faces already introduced above    !
    !---------------------------------------------------------!
    do e = 1, Prim % n_edges
      if(edge_data(e) .ne. 0) then

        ! Take the Prim edge's nodes (these are from Prim)
        do i_nod = 1, 2
          n_p = Prim % edges_n(i_nod, e)

          ! This node_to_face was stored in the previous step
          f_d = node_to_face(n_p)
          Dual % faces_n_nodes(f_d) = Dual % faces_n_nodes(f_d) + 1
          call Adjust_First_Dim(Dual % faces_n_nodes(f_d), Dual % faces_n)
          Dual % faces_n(Dual % faces_n_nodes(f_d), f_d) = edge_to_node(e)

          ! This node_to_cell was stored in the previous step
          b_d = node_to_cell(n_p)
          Dual % cells_n_nodes(b_d) = Dual % cells_n_nodes(b_d) + 1
          call Adjust_First_Dim(Dual % cells_n_nodes(b_d), Dual % cells_n)
          Dual % cells_n(Dual % cells_n_nodes(b_d), b_d) = edge_to_node(e)
        end do  ! i_nod for edge, goes from 1 to 2

      end if
    end do  ! through edges

    !------------------------------------------------------------!
    !   Add additional nodes to Dual's face from sharp corners   !
    !    Here we work on the faces already introduced before     !
    !------------------------------------------------------------!
    do e = 1, Prim % n_edges
      if(edge_data(e) .ne. 0) then

        do i_nod = 1, 2
          n_p = Prim % edges_n(i_nod, e)
          f_d = node_to_face(n_p)
          b_d = node_to_cell(n_p)

          ! The grid has sharp corners, add them to boundary faces and cells
          if(sharp_corner(n_p) .gt. 0) then
            if(sharp_inject(f_d) .eq. 0) then

              ! Estimate the number of new node in Dual
              ! (sharp_corner holds its local number)
              n_d = Dual % n_nodes - sharp_corner(n_p) + 1
              Dual % xn(n_d) = Prim % xn(n_p)
              Dual % yn(n_d) = Prim % yn(n_p)
              Dual % zn(n_d) = Prim % zn(n_p)

              Dual % faces_n_nodes(f_d) = Dual % faces_n_nodes(f_d) + 1
              call Adjust_First_Dim(Dual % faces_n_nodes(f_d), Dual % faces_n)
              Dual % faces_n(Dual % faces_n_nodes(f_d), f_d) = n_d

              Dual % cells_n_nodes(b_d) = Dual % cells_n_nodes(b_d) + 1
              call Adjust_First_Dim(Dual % cells_n_nodes(b_d), Dual % cells_n)
              Dual % cells_n(Dual % cells_n_nodes(b_d), b_d) = n_d

              ! Mark that the face has been injected a sharp corner
              ! write(*,'(a,i9)') ' # Injecting new node in face', f_d
              sharp_inject(f_d) = sharp_inject(f_d) + 1

              ! Check those little nodes inserted just before
              n1 = Dual % faces_n(Dual % faces_n_nodes(f_d)-2, f_d)
              n2 = Dual % faces_n(Dual % faces_n_nodes(f_d)-1, f_d)
              if(concave_link(1,n_d) .eq. 0 .and.  &
                 concave_link(2,n_d) .eq. 0) then
                concave_link(1,n_d) = n1
                concave_link(2,n_d) = n2
              else
                print *, '# Well, you didn''t see this coming!'
                print *, '# One sharp corner can be in more than one edge'
              end if

            end if  ! sharp inject in f_d
          end if    ! sharp corner here

        end do  ! i_nod for edge, goes from 1 to 2
      end if
    end do  ! through edges

  end do

  !---------------------------------------!
  !                                       !
  !   Sort the points on boundary faces   !
  !                                       !
  !---------------------------------------!
  do s = Prim % n_edges + 1, Dual % n_faces
    call Convert % Sort_Face_Nodes(Dual, s, concave_link)
  end do

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
  do s = 1, Dual % n_faces
    c1 = Dual % faces_c(1, s)
    c2 = Dual % faces_c(2, s)

    !----------------------------------------------------!
    !   Store faces surrounding each cell in Dual grid   !
    !----------------------------------------------------!
    Dual % cells_n_faces(c1) = Dual % cells_n_faces(c1) + 1
    Dual % cells_n_faces(c2) = Dual % cells_n_faces(c2) + 1
    call Adjust_First_Dim(Dual % cells_n_faces(c1), Dual % cells_f)
    call Adjust_First_Dim(Dual % cells_n_faces(c2), Dual % cells_f)
    Dual % cells_f(Dual % cells_n_faces(c1), c1) = s
    Dual % cells_f(Dual % cells_n_faces(c2), c2) = s

    !----------------------------------------------------!
    !   Store nodes surrounding each cell in Dual grid   !
    !----------------------------------------------------!
    do i_nod = 1, Dual % faces_n_nodes(s)
      n = Dual % faces_n(i_nod, s)

      ! Handle cell 1
      do j_nod = 1, Dual % cells_n_nodes(c1)
        n1 = Dual % cells_n(j_nod, c1)
        if(n1 .eq. n) goto 1
      end do  ! j_nod
      Dual % cells_n_nodes(c1) = Dual % cells_n_nodes(c1) + 1
      call Adjust_First_Dim(Dual % cells_n_nodes(c1), Dual % cells_n)
      Dual % cells_n(Dual % cells_n_nodes(c1), c1) = n
1     continue

      ! Handle cell 2
      do j_nod = 1, Dual % cells_n_nodes(c2)
        n2 = Dual % cells_n(j_nod, c2)
        if(n2 .eq. n) goto 2
      end do  ! j_nod
      Dual % cells_n_nodes(c2) = Dual % cells_n_nodes(c2) + 1
      call Adjust_First_Dim(Dual % cells_n_nodes(c2), Dual % cells_n)
      Dual % cells_n(Dual % cells_n_nodes(c2), c2) = n
2     continue

    end do  ! do i_nod

  end do  ! do s

  !--------------------------------------------------------!
  !                                                        !
  !   Fix the number of nodes for polyhedral (all) cells   !
  !                                                        !
  !--------------------------------------------------------!
  do c = 1, Dual % n_cells
    call Sort % Int_Array(Dual % cells_n(1:Dual % cells_n_nodes(c),c))
    Dual % cells_n_nodes(c) = -Dual % cells_n_nodes(c)
  end do

  end subroutine
