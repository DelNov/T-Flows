!==============================================================================!
  subroutine Refines_Mod_Refine_Marked_Cells(ref, Grid, lev)
!------------------------------------------------------------------------------!
!>  This function is responsible for refining the grid cells that have been
!>  marked for refinement. It adjusts the grid structure to account for the
!>  increased resolution in specific areas, adding new nodes and cells as
!>  required.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Cell and node handling: It processes each cell in the grid. For the      !
!     cells marked for refinement, it subdivides each cell into eight smaller  !
!     cells. This involves calculating the positions of new nodes (such as     !
!     edge midpoints, face centers, and cell center) and updating the          !
!     cell-to-node connectivity.                                               !
!   * Neighbor update: For each refined cell, it updates the neighboring cell  !
!     information, ensuring the grid's integrity and connectivity.             !
!   * Refinement levels: The refinement level for each new cell is set,        !
!     keeping track of the depth of refinement.                                !
!   * Handling new nodes: The subroutine calculates the positions of new nodes !
!     by averaging the positions of existing nodes. These new nodes include    !
!     those placed on the edges, faces, and center of the original cell.       !
!   * Connectivity update: It updates the connectivity information for cells   !
!     and nodes, including handling the special cases where neighboring cells  !
!     are also refined.                                                        !
!   * Cleanup and rearrangement: After refinement, it cleans up and rearranges !
!     the data structures to reflect the new grid configuration. This includes !
!     renumbering cells and removing redundant data.                           !
!   * Memory management: Allocation and deallocation of temporary arrays used  !
!     for storing node and cell information during the refinement process.     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Refines_Type) :: ref
  type(Grid_Type)    :: Grid
  integer            :: lev
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, n_cells_old, c1, c2, c3, c4, c5, c6
  integer :: cr1, cr2, cr3, cr4, cr5, cr6, cr7, cr8
  integer :: n, n_nodes_old, n1, n2, n3, n4, n5, n6, n7, n8
  integer :: n12, n13, n24, n34, n15, n26, n37, n48, n56, n57, n68, n78
  integer :: nf1, nf2, nf3, nf4, nf5, nf6, n0
  integer :: del   ! number of deleted nodes
  integer :: na, na0, na1, na2, nb, nb0, nb1, nb2
  integer :: nn_2, nn_4, nn_8
  integer, allocatable :: node_n2(:,:)
  integer, allocatable :: node_n4(:,:)
  integer, allocatable :: node_n8(:,:)
  integer, allocatable :: cell_points_to(:)
!==============================================================================!
!                                                                              !
!                               c6      c3                                     !
!                               |      /                                       !
!                         8-----|---------6                                    !
!                        /  cr8 |/  cr6  /|                                    !
!                       /-------+-------/ |                                    !
!                      /       /       /| |                                    !
!                     7-------+-------5 | |                                    !
!                     |       |       | |/|                                    !
!                c4---|  cr7  |  cr5  | +-----c2                               !
!                     |       |       |/| |                                    !
!                     +-------+-------+ | 2                                    !
!                     |      /|       | |/                                     !
!                     |  cr3/ |  cr1  | /                                      !
!                     |    /  |       |/                                       !
!                     3---c5--+-------1                                        !
!                               |                                              !
!                               c1                                             !
!                                                                              !
!------------------------------------------------------------------------------!

  print *, '#=============================================='
  print *, '# Refine: Number of nodes: ', Grid % n_nodes
  print *, '#         Number of cells: ', Grid % n_cells
  print *, '#----------------------------------------------'

  n_cells_old = Grid % n_cells
  n_nodes_old = Grid % n_nodes
  nn_2 = 0
  nn_4 = 0
  nn_8 = 0

  ! Allocate local memory
  call Enlarge % Matrix_Int(node_n2, i=(/1,Grid % n_nodes/), j=(/0,2/))
  call Enlarge % Matrix_Int(node_n4, i=(/1,Grid % n_nodes/), j=(/0,4/))
  call Enlarge % Matrix_Int(node_n8, i=(/1,Grid % n_nodes/), j=(/0,8/))
  call Enlarge % Array_Int(cell_points_to, i=(/1,Grid % n_nodes/))

  !---------------------!
  !                     !
  !   Count new celss   !
  !                     !
  !---------------------!
  do c = 1, n_cells_old
    if(ref % cell_marked(c)) then
      Grid % n_cells = Grid % n_cells + 8
      cell_points_to(c) = Grid % n_cells   ! now points to cr8
    end if
  end do

  ! Expand cell arrays to accomodate new cells
  ! (At this stage, only Grid % n_cells matters)
  call Grid % Allocate_Cells(Grid % n_cells, Grid % n_bnd_cells)
  call Enlarge % Array_Int(cell_points_to, i=(/1,Grid % n_nodes/))
  call Refines_Mod_Allocate_Cells(ref, Grid % n_cells, Grid % n_bnd_cells)

  do c = 1, n_cells_old

    c1 = Grid % cells_c(1,c)
    c2 = Grid % cells_c(2,c)
    c3 = Grid % cells_c(3,c)
    c4 = Grid % cells_c(4,c)
    c5 = Grid % cells_c(5,c)
    c6 = Grid % cells_c(6,c)

    n1 = Grid % cells_n(1,c)
    n2 = Grid % cells_n(2,c)
    n3 = Grid % cells_n(3,c)
    n4 = Grid % cells_n(4,c)
    n5 = Grid % cells_n(5,c)
    n6 = Grid % cells_n(6,c)
    n7 = Grid % cells_n(7,c)
    n8 = Grid % cells_n(8,c)

    !-------------------!
    !                   !
    !   Refined cells   !
    !                   !
    !-------------------!
    if( ref % cell_marked(c) ) then  ! only refined; ends at line 826

      !------------------------------------!
      !   Take care of neighboring cells   !
      !------------------------------------!
      cr1 = cell_points_to(c) - 7
      cr2 = cell_points_to(c) - 6
      cr3 = cell_points_to(c) - 5
      cr4 = cell_points_to(c) - 4
      cr5 = cell_points_to(c) - 3
      cr6 = cell_points_to(c) - 2
      cr7 = cell_points_to(c) - 1
      cr8 = cell_points_to(c)

      !-----------------------------------------------!
      !   Internal links do not depend on neighbors   !
      !-----------------------------------------------!

      ! 6
      Grid % cells_c(6,cr1) = cr5
      Grid % cells_c(6,cr2) = cr6
      Grid % cells_c(6,cr3) = cr7
      Grid % cells_c(6,cr4) = cr8

      ! 5
      Grid % cells_c(5,cr2) = cr1
      Grid % cells_c(5,cr4) = cr3
      Grid % cells_c(5,cr6) = cr5
      Grid % cells_c(5,cr8) = cr7

      ! 4
      Grid % cells_c(4,cr1) = cr3
      Grid % cells_c(4,cr2) = cr4
      Grid % cells_c(4,cr5) = cr7
      Grid % cells_c(4,cr6) = cr8

      ! 3
      Grid % cells_c(3,cr1) = cr2
      Grid % cells_c(3,cr3) = cr4
      Grid % cells_c(3,cr5) = cr6
      Grid % cells_c(3,cr7) = cr8

      ! 2
      Grid % cells_c(2,cr3) = cr1
      Grid % cells_c(2,cr4) = cr2
      Grid % cells_c(2,cr7) = cr5
      Grid % cells_c(2,cr8) = cr6

      ! 1
      Grid % cells_c(1,cr5) = cr1
      Grid % cells_c(1,cr6) = cr2
      Grid % cells_c(1,cr7) = cr3
      Grid % cells_c(1,cr8) = cr4

      !-------------------------!
      !   Level of refinement   !
      !-------------------------!
      ref % cell_level(cr1) = lev
      ref % cell_level(cr2) = lev
      ref % cell_level(cr3) = lev
      ref % cell_level(cr4) = lev
      ref % cell_level(cr5) = lev
      ref % cell_level(cr6) = lev
      ref % cell_level(cr7) = lev
      ref % cell_level(cr8) = lev

      !----------------------------------------!
      !   External links depend on neighbors   !
      !----------------------------------------!
      if(.not. ref % cell_marked(c1)) then  ! neighbor 1 not refined
        Grid % cells_c(1,cr1) = c1
        Grid % cells_c(1,cr2) = c1
        Grid % cells_c(1,cr3) = c1
        Grid % cells_c(1,cr4) = c1
      else                        ! neighbor 1 refined
        Grid % cells_c(1,cr1) = cell_points_to(c1)-8   &
                              + Refines_Mod_Which_Node(ref,Grid,c1,n1)
        Grid % cells_c(1,cr2) = cell_points_to(c1)-8   &
                              + Refines_Mod_Which_Node(ref,Grid,c1,n2)
        Grid % cells_c(1,cr3) = cell_points_to(c1)-8   &
                              + Refines_Mod_Which_Node(ref,Grid,c1,n3)
        Grid % cells_c(1,cr4) = cell_points_to(c1)-8   &
                              + Refines_Mod_Which_Node(ref,Grid,c1,n4)
      end if

      if(.not. ref % cell_marked(c2)) then  ! neighbor 2 not refined
        Grid % cells_c(2,cr1) = c2
        Grid % cells_c(2,cr2) = c2
        Grid % cells_c(2,cr5) = c2
        Grid % cells_c(2,cr6) = c2
      else                        ! neighbor 2 refined
        Grid % cells_c(2,cr1) = cell_points_to(c2)-8   &
                              + Refines_Mod_Which_Node(ref,Grid,c2,n1)
        Grid % cells_c(2,cr2) = cell_points_to(c2)-8   &
                              + Refines_Mod_Which_Node(ref,Grid,c2,n2)
        Grid % cells_c(2,cr5) = cell_points_to(c2)-8   &
                              + Refines_Mod_Which_Node(ref,Grid,c2,n5)
        Grid % cells_c(2,cr6) = cell_points_to(c2)-8   &
                              + Refines_Mod_Which_Node(ref,Grid,c2,n6)
      end if

      if(.not. ref % cell_marked(c3)) then  ! neighbor 3 not refined
        Grid % cells_c(3,cr2) = c3
        Grid % cells_c(3,cr4) = c3
        Grid % cells_c(3,cr6) = c3
        Grid % cells_c(3,cr8) = c3
      else                        ! neighbor 3 refined
        Grid % cells_c(3,cr2) = cell_points_to(c3)-8   &
                              + Refines_Mod_Which_Node(ref,Grid,c3,n2)
        Grid % cells_c(3,cr4) = cell_points_to(c3)-8   &
                              + Refines_Mod_Which_Node(ref,Grid,c3,n4)
        Grid % cells_c(3,cr6) = cell_points_to(c3)-8   &
                              + Refines_Mod_Which_Node(ref,Grid,c3,n6)
        Grid % cells_c(3,cr8) = cell_points_to(c3)-8   &
                              + Refines_Mod_Which_Node(ref,Grid,c3,n8)
      end if

      if(.not. ref % cell_marked(c4)) then  ! neighbor 4 not refine
        Grid % cells_c(4,cr3) = c4
        Grid % cells_c(4,cr4) = c4
        Grid % cells_c(4,cr7) = c4
        Grid % cells_c(4,cr8) = c4
      else                        ! neighbor 4 refine
        Grid % cells_c(4,cr3) = cell_points_to(c4)-8   &
                              + Refines_Mod_Which_Node(ref,Grid,c4,n3)
        Grid % cells_c(4,cr4) = cell_points_to(c4)-8   &
                              + Refines_Mod_Which_Node(ref,Grid,c4,n4)
        Grid % cells_c(4,cr7) = cell_points_to(c4)-8   &
                              + Refines_Mod_Which_Node(ref,Grid,c4,n7)
        Grid % cells_c(4,cr8) = cell_points_to(c4)-8   &
                              + Refines_Mod_Which_Node(ref,Grid,c4,n8)
      end if

      if(.not. ref % cell_marked(c5)) then  ! neighbor 5 not refined
        Grid % cells_c(5,cr1) = c5
        Grid % cells_c(5,cr3) = c5
        Grid % cells_c(5,cr5) = c5
        Grid % cells_c(5,cr7) = c5
      else                        ! neighbor 5 refined
        Grid % cells_c(5,cr1) = cell_points_to(c5)-8   &
                              + Refines_Mod_Which_Node(ref,Grid,c5,n1)
        Grid % cells_c(5,cr3) = cell_points_to(c5)-8   &
                              + Refines_Mod_Which_Node(ref,Grid,c5,n3)
        Grid % cells_c(5,cr5) = cell_points_to(c5)-8   &
                              + Refines_Mod_Which_Node(ref,Grid,c5,n5)
        Grid % cells_c(5,cr7) = cell_points_to(c5)-8   &
                              + Refines_Mod_Which_Node(ref,Grid,c5,n7)
      end if

      if(.not. ref % cell_marked(c6)) then  ! neighbor 6 not refined
        Grid % cells_c(6,cr5) = c6
        Grid % cells_c(6,cr6) = c6
        Grid % cells_c(6,cr7) = c6
        Grid % cells_c(6,cr8) = c6
      else                        ! neighbor 6 refined
        Grid % cells_c(6,cr5) = cell_points_to(c6)-8   &
                              + Refines_Mod_Which_Node(ref,Grid,c6,n5)
        Grid % cells_c(6,cr6) = cell_points_to(c6)-8   &
                              + Refines_Mod_Which_Node(ref,Grid,c6,n6)
        Grid % cells_c(6,cr7) = cell_points_to(c6)-8   &
                              + Refines_Mod_Which_Node(ref,Grid,c6,n7)
        Grid % cells_c(6,cr8) = cell_points_to(c6)-8   &
                              + Refines_Mod_Which_Node(ref,Grid,c6,n8)
      end if

      !------------------------!
      !   Take care of nodes   !
      !------------------------!

      !------------------------!
      !   Nodes on the edges   !
      !------------------------!
      n12 = 0     !
      n13 = 0     !         8-----n68-----6
      n24 = 0     !        /|            /|
      n34 = 0     !      n78|          n56|
      n15 = 0     !      / n48         / n26
      n26 = 0     !     7-----n57-----5   |
      n37 = 0     !     |   |         |   |
      n48 = 0     !     |   4- - -n24-| - 2
      n56 = 0     !    n37 /         n15 /
      n57 = 0     !     |n34          |n12
      n68 = 0     !     |/            |/
      n78 = 0     !     3-----n13-----1

      ! n12
      do n = 1, nn_2
        if( ( node_n2(n,1) .eq. n1 .and. node_n2(n,2) .eq. n2 ) .or.  &
            ( node_n2(n,1) .eq. n2 .and. node_n2(n,2) .eq. n1) ) then
          n12 = node_n2(n,0)
        end if
      end do
      if (n12 .eq. 0) then
        nn_2 = nn_2 + 1
        Grid % n_nodes = Grid % n_nodes  + 1
        call Grid % Allocate_Nodes(Grid % n_nodes, GROWTH_MARGIN)
        n12 = Grid % n_nodes
        node_n2(nn_2,0) = n12
        node_n2(nn_2,1) = n1
        node_n2(nn_2,2) = n2
        Grid % xn(n12) = .5 * (Grid % xn(n1) + Grid % xn(n2))
        Grid % yn(n12) = .5 * (Grid % yn(n1) + Grid % yn(n2))
        Grid % zn(n12) = .5 * (Grid % zn(n1) + Grid % zn(n2))
      end if

      ! n13
      do n = 1, nn_2
        if( ( node_n2(n,1) .eq. n1 .and. node_n2(n,2) .eq. n3 ) .or.  &
            ( node_n2(n,1) .eq. n3 .and. node_n2(n,2) .eq. n1) ) then
          n13 = node_n2(n,0)
        end if
      end do
      if (n13 .eq. 0) then
        nn_2 = nn_2 + 1
        Grid % n_nodes = Grid % n_nodes  + 1
        call Grid % Allocate_Nodes(Grid % n_nodes, GROWTH_MARGIN)
        n13 = Grid % n_nodes
        node_n2(nn_2,0) = n13
        node_n2(nn_2,1) = n1
        node_n2(nn_2,2) = n3
        Grid % xn(n13) = .5 * (Grid % xn(n1) + Grid % xn(n3))
        Grid % yn(n13) = .5 * (Grid % yn(n1) + Grid % yn(n3))
        Grid % zn(n13) = .5 * (Grid % zn(n1) + Grid % zn(n3))
      end if

      ! n24
      do n = 1, nn_2
        if( ( node_n2(n,1) .eq. n2 .and. node_n2(n,2) .eq. n4 ) .or.  &
            ( node_n2(n,1) .eq. n4 .and. node_n2(n,2) .eq. n2) ) then
          n24 = node_n2(n,0)
        end if
      end do
      if (n24 .eq. 0) then
        nn_2 = nn_2 + 1
        Grid % n_nodes = Grid % n_nodes  + 1
        call Grid % Allocate_Nodes(Grid % n_nodes, GROWTH_MARGIN)
        n24 = Grid % n_nodes
        node_n2(nn_2,0) = n24
        node_n2(nn_2,1) = n2
        node_n2(nn_2,2) = n4
        Grid % xn(n24) = .5 * (Grid % xn(n2) + Grid % xn(n4))
        Grid % yn(n24) = .5 * (Grid % yn(n2) + Grid % yn(n4))
        Grid % zn(n24) = .5 * (Grid % zn(n2) + Grid % zn(n4))
      end if

      ! n34
      do n = 1, nn_2
        if( ( node_n2(n,1) .eq. n3 .and. node_n2(n,2) .eq. n4 ) .or.  &
            ( node_n2(n,1) .eq. n4 .and. node_n2(n,2) .eq. n3) ) then
          n34 = node_n2(n,0)
        end if
      end do
      if (n34 .eq. 0) then
        nn_2 = nn_2 + 1
        Grid % n_nodes = Grid % n_nodes  + 1
        call Grid % Allocate_Nodes(Grid % n_nodes, GROWTH_MARGIN)
        n34 = Grid % n_nodes
        node_n2(nn_2,0) = n34
        node_n2(nn_2,1) = n3
        node_n2(nn_2,2) = n4
        Grid % xn(n34) = .5 * (Grid % xn(n3) + Grid % xn(n4))
        Grid % yn(n34) = .5 * (Grid % yn(n3) + Grid % yn(n4))
        Grid % zn(n34) = .5 * (Grid % zn(n3) + Grid % zn(n4))
      end if

      ! n15
      do n = 1, nn_2
        if( ( node_n2(n,1) .eq. n1 .and. node_n2(n,2) .eq. n5 ) .or.  &
            ( node_n2(n,1) .eq. n5 .and. node_n2(n,2) .eq. n1) ) then
          n15 = node_n2(n,0)
        end if
      end do
      if (n15 .eq. 0) then
        nn_2 = nn_2 + 1
        Grid % n_nodes = Grid % n_nodes  + 1
        call Grid % Allocate_Nodes(Grid % n_nodes, GROWTH_MARGIN)
        n15 = Grid % n_nodes
        node_n2(nn_2,0) = n15
        node_n2(nn_2,1) = n1
        node_n2(nn_2,2) = n5
        Grid % xn(n15) = .5 * (Grid % xn(n1) + Grid % xn(n5))
        Grid % yn(n15) = .5 * (Grid % yn(n1) + Grid % yn(n5))
        Grid % zn(n15) = .5 * (Grid % zn(n1) + Grid % zn(n5))
      end if

      ! n26
      do n = 1, nn_2
        if( ( node_n2(n,1) .eq. n2 .and. node_n2(n,2) .eq. n6 ) .or.  &
            ( node_n2(n,1) .eq. n6 .and. node_n2(n,2) .eq. n2) ) then
          n26 = node_n2(n,0)
        end if
      end do
      if (n26 .eq. 0) then
        nn_2 = nn_2 + 1
        Grid % n_nodes = Grid % n_nodes  + 1
        call Grid % Allocate_Nodes(Grid % n_nodes, GROWTH_MARGIN)
        n26 = Grid % n_nodes
        node_n2(nn_2,0) = n26
        node_n2(nn_2,1) = n2
        node_n2(nn_2,2) = n6
        Grid % xn(n26) = .5 * (Grid % xn(n2) + Grid % xn(n6))
        Grid % yn(n26) = .5 * (Grid % yn(n2) + Grid % yn(n6))
        Grid % zn(n26) = .5 * (Grid % zn(n2) + Grid % zn(n6))
      end if

      ! n37
      do n = 1, nn_2
        if( ( node_n2(n,1) .eq. n3 .and. node_n2(n,2) .eq. n7 ) .or.  &
            ( node_n2(n,1) .eq. n7 .and. node_n2(n,2) .eq. n3) ) then
          n37 = node_n2(n,0)
        end if
      end do
      if (n37 .eq. 0) then
        nn_2 = nn_2 + 1
        Grid % n_nodes = Grid % n_nodes  + 1
        call Grid % Allocate_Nodes(Grid % n_nodes, GROWTH_MARGIN)
        n37 = Grid % n_nodes
        node_n2(nn_2,0) = n37
        node_n2(nn_2,1) = n3
        node_n2(nn_2,2) = n7
        Grid % xn(n37) = .5 * (Grid % xn(n3) + Grid % xn(n7))
        Grid % yn(n37) = .5 * (Grid % yn(n3) + Grid % yn(n7))
        Grid % zn(n37) = .5 * (Grid % zn(n3) + Grid % zn(n7))
      end if

      ! n48
      do n = 1, nn_2
        if( ( node_n2(n,1) .eq. n4 .and. node_n2(n,2) .eq. n8 ) .or.  &
            ( node_n2(n,1) .eq. n8 .and. node_n2(n,2) .eq. n4) ) then
          n48 = node_n2(n,0)
        end if
      end do
      if (n48 .eq. 0) then
        nn_2 = nn_2 + 1
        Grid % n_nodes = Grid % n_nodes  + 1
        call Grid % Allocate_Nodes(Grid % n_nodes, GROWTH_MARGIN)
        n48 = Grid % n_nodes
        node_n2(nn_2,0) = n48
        node_n2(nn_2,1) = n4
        node_n2(nn_2,2) = n8
        Grid % xn(n48) = .5 * (Grid % xn(n4) + Grid % xn(n8))
        Grid % yn(n48) = .5 * (Grid % yn(n4) + Grid % yn(n8))
        Grid % zn(n48) = .5 * (Grid % zn(n4) + Grid % zn(n8))
      end if

      ! n56
      do n = 1, nn_2
        if( ( node_n2(n,1) .eq. n5 .and. node_n2(n,2) .eq. n6 ) .or.  &
            ( node_n2(n,1) .eq. n6 .and. node_n2(n,2) .eq. n5) ) then
          n56 = node_n2(n,0)
        end if
      end do
      if (n56 .eq. 0) then
        nn_2 = nn_2 + 1
        Grid % n_nodes = Grid % n_nodes  + 1
        call Grid % Allocate_Nodes(Grid % n_nodes, GROWTH_MARGIN)
        n56 = Grid % n_nodes
        node_n2(nn_2,0) = n56
        node_n2(nn_2,1) = n5
        node_n2(nn_2,2) = n6
        Grid % xn(n56) = .5 * (Grid % xn(n5) + Grid % xn(n6))
        Grid % yn(n56) = .5 * (Grid % yn(n5) + Grid % yn(n6))
        Grid % zn(n56) = .5 * (Grid % zn(n5) + Grid % zn(n6))
      end if

      ! n57
      do n = 1, nn_2
        if( ( node_n2(n,1) .eq. n5 .and. node_n2(n,2) .eq. n7 ) .or.  &
            ( node_n2(n,1) .eq. n7 .and. node_n2(n,2) .eq. n5) ) then
          n57 = node_n2(n,0)
        end if
      end do
      if (n57 .eq. 0) then
        nn_2 = nn_2 + 1
        Grid % n_nodes = Grid % n_nodes  + 1
        call Grid % Allocate_Nodes(Grid % n_nodes, GROWTH_MARGIN)
        n57 = Grid % n_nodes
        node_n2(nn_2,0) = n57
        node_n2(nn_2,1) = n5
        node_n2(nn_2,2) = n7
        Grid % xn(n57) = .5 * (Grid % xn(n5) + Grid % xn(n7))
        Grid % yn(n57) = .5 * (Grid % yn(n5) + Grid % yn(n7))
        Grid % zn(n57) = .5 * (Grid % zn(n5) + Grid % zn(n7))
      end if

      ! n68
      do n = 1, nn_2
        if( ( node_n2(n,1) .eq. n6 .and. node_n2(n,2) .eq. n8 ) .or.  &
            ( node_n2(n,1) .eq. n8 .and. node_n2(n,2) .eq. n6) ) then
          n68 = node_n2(n,0)
        end if
      end do
      if (n68 .eq. 0) then
        nn_2 = nn_2 + 1
        Grid % n_nodes = Grid % n_nodes  + 1
        call Grid % Allocate_Nodes(Grid % n_nodes, GROWTH_MARGIN)
        n68 = Grid % n_nodes
        node_n2(nn_2,0) = n68
        node_n2(nn_2,1) = n6
        node_n2(nn_2,2) = n8
        Grid % xn(n68) = .5 * (Grid % xn(n6) + Grid % xn(n8))
        Grid % yn(n68) = .5 * (Grid % yn(n6) + Grid % yn(n8))
        Grid % zn(n68) = .5 * (Grid % zn(n6) + Grid % zn(n8))
      end if

      ! n78
      do n = 1, nn_2
        if( ( node_n2(n,1) .eq. n7 .and. node_n2(n,2) .eq. n8 ) .or.  &
            ( node_n2(n,1) .eq. n8 .and. node_n2(n,2) .eq. n7) ) then
          n78 = node_n2(n,0)
        end if
      end do
      if (n78 .eq. 0) then
        nn_2 = nn_2 + 1
        Grid % n_nodes = Grid % n_nodes  + 1
        call Grid % Allocate_Nodes(Grid % n_nodes, GROWTH_MARGIN)
        n78 = Grid % n_nodes
        node_n2(nn_2,0) = n78
        node_n2(nn_2,1) = n7
        node_n2(nn_2,2) = n8
        Grid % xn(n78) = .5 * (Grid % xn(n7) + Grid % xn(n8))
        Grid % yn(n78) = .5 * (Grid % yn(n7) + Grid % yn(n8))
        Grid % zn(n78) = .5 * (Grid % zn(n7) + Grid % zn(n8))
      end if

      !-------------------------!
      !   Then nodes on faces   !
      !-------------------------!
      nf1 = 0
      nf2 = 0
      nf3 = 0
      nf4 = 0
      nf5 = 0
      nf6 = 0

      ! nf1
      do n = 1, nn_4
        if( ( node_n4(n,1) .eq. n1 .and. node_n4(n,4) .eq. n4 ) .or.  &
            ( node_n4(n,1) .eq. n4 .and. node_n4(n,4) .eq. n1 ) .or.  &
            ( node_n4(n,1) .eq. n2 .and. node_n4(n,4) .eq. n3 ) .or.  &
            ( node_n4(n,1) .eq. n3 .and. node_n4(n,4) .eq. n2 ) ) then
          nf1 = node_n4(n,0)
        end if
      end do
      if (nf1 .eq. 0) then
        nn_4 = nn_4 + 1
        Grid % n_nodes = Grid % n_nodes  + 1
        call Grid % Allocate_Nodes(Grid % n_nodes, GROWTH_MARGIN)
        nf1 = Grid % n_nodes
        node_n4(nn_4,0) = nf1
        node_n4(nn_4,1) = n1
        node_n4(nn_4,2) = n2
        node_n4(nn_4,3) = n3
        node_n4(nn_4,4) = n4
        Grid % xn(nf1) = 0.25 * (Grid % xn(n1) + Grid % xn(n2) +  &
                                 Grid % xn(n3) + Grid % xn(n4))
        Grid % yn(nf1) = 0.25 * (Grid % yn(n1) + Grid % yn(n2) +  &
                                 Grid % yn(n3) + Grid % yn(n4))
        Grid % zn(nf1) = 0.25 * (Grid % zn(n1) + Grid % zn(n2) +  &
                                 Grid % zn(n3) + Grid % zn(n4))
      end if

      ! nf2
      do n = 1, nn_4
        if( ( node_n4(n,1) .eq. n1 .and. node_n4(n,4) .eq. n6 ) .or.  &
            ( node_n4(n,1) .eq. n6 .and. node_n4(n,4) .eq. n1 ) .or.  &
            ( node_n4(n,1) .eq. n2 .and. node_n4(n,4) .eq. n5 ) .or.  &
            ( node_n4(n,1) .eq. n5 .and. node_n4(n,4) .eq. n2 ) ) then
          nf2 = node_n4(n,0)
        end if
      end do
      if (nf2 .eq. 0) then
        nn_4 = nn_4 + 1
        Grid % n_nodes = Grid % n_nodes  + 1
        call Grid % Allocate_Nodes(Grid % n_nodes, GROWTH_MARGIN)
        nf2 = Grid % n_nodes
        node_n4(nn_4,0) = nf2
        node_n4(nn_4,1) = n1
        node_n4(nn_4,2) = n2
        node_n4(nn_4,3) = n5
        node_n4(nn_4,4) = n6
        Grid % xn(nf2) = 0.25 * (Grid % xn(n1) + Grid % xn(n2) +  &
                                 Grid % xn(n5) + Grid % xn(n6))
        Grid % yn(nf2) = 0.25 * (Grid % yn(n1) + Grid % yn(n2) +  &
                                 Grid % yn(n5) + Grid % yn(n6))
        Grid % zn(nf2) = 0.25 * (Grid % zn(n1) + Grid % zn(n2) +  &
                                 Grid % zn(n5) + Grid % zn(n6))
      end if

      ! nf3
      do n = 1, nn_4
        if( ( node_n4(n,1) .eq. n2 .and. node_n4(n,4) .eq. n8 ) .or.  &
            ( node_n4(n,1) .eq. n8 .and. node_n4(n,4) .eq. n2 ) .or.  &
            ( node_n4(n,1) .eq. n4 .and. node_n4(n,4) .eq. n6 ) .or.  &
            ( node_n4(n,1) .eq. n6 .and. node_n4(n,4) .eq. n4 ) ) then
          nf3 = node_n4(n,0)
        end if
      end do
      if (nf3 .eq. 0) then
        nn_4 = nn_4 + 1
        Grid % n_nodes = Grid % n_nodes  + 1
        call Grid % Allocate_Nodes(Grid % n_nodes, GROWTH_MARGIN)
        nf3 = Grid % n_nodes
        node_n4(nn_4,0) = nf3
        node_n4(nn_4,1) = n2
        node_n4(nn_4,2) = n4
        node_n4(nn_4,3) = n6
        node_n4(nn_4,4) = n8
        Grid % xn(nf3) = 0.25 * (Grid % xn(n2) + Grid % xn(n4) +  &
                                 Grid % xn(n6) + Grid % xn(n8))
        Grid % yn(nf3) = 0.25 * (Grid % yn(n2) + Grid % yn(n4) +  &
                                 Grid % yn(n6) + Grid % yn(n8))
        Grid % zn(nf3) = 0.25 * (Grid % zn(n2) + Grid % zn(n4) +  &
                                 Grid % zn(n6) + Grid % zn(n8))
      end if

      ! nf4
      do n = 1, nn_4
        if( ( node_n4(n,1) .eq. n3 .and. node_n4(n,4) .eq. n8 ) .or.  &
            ( node_n4(n,1) .eq. n8 .and. node_n4(n,4) .eq. n3 ) .or.  &
            ( node_n4(n,1) .eq. n4 .and. node_n4(n,4) .eq. n7 ) .or.  &
            ( node_n4(n,1) .eq. n7 .and. node_n4(n,4) .eq. n4 ) ) then
          nf4 = node_n4(n,0)
        end if
      end do
      if (nf4 .eq. 0) then
        nn_4 = nn_4 + 1
        Grid % n_nodes = Grid % n_nodes  + 1
        call Grid % Allocate_Nodes(Grid % n_nodes, GROWTH_MARGIN)
        nf4 = Grid % n_nodes
        node_n4(nn_4,0) = nf4
        node_n4(nn_4,1) = n3
        node_n4(nn_4,2) = n4
        node_n4(nn_4,3) = n7
        node_n4(nn_4,4) = n8
        Grid % xn(nf4) = 0.25 * (Grid % xn(n3) + Grid % xn(n4) +  &
                                 Grid % xn(n7) + Grid % xn(n8))
        Grid % yn(nf4) = 0.25 * (Grid % yn(n3) + Grid % yn(n4) +  &
                                 Grid % yn(n7) + Grid % yn(n8))
        Grid % zn(nf4) = 0.25 * (Grid % zn(n3) + Grid % zn(n4) +  &
                                 Grid % zn(n7) + Grid % zn(n8))
      end if

      ! nf5
      do n = 1, nn_4
        if( ( node_n4(n,1) .eq. n1 .and. node_n4(n,4) .eq. n7 ) .or.  &
            ( node_n4(n,1) .eq. n7 .and. node_n4(n,4) .eq. n1 ) .or.  &
            ( node_n4(n,1) .eq. n3 .and. node_n4(n,4) .eq. n5 ) .or.  &
            ( node_n4(n,1) .eq. n5 .and. node_n4(n,4) .eq. n3 ) ) then
          nf5 = node_n4(n,0)
        end if
      end do
      if (nf5 .eq. 0) then
        nn_4 = nn_4 + 1
        Grid % n_nodes = Grid % n_nodes  + 1
        call Grid % Allocate_Nodes(Grid % n_nodes, GROWTH_MARGIN)
        nf5 = Grid % n_nodes
        node_n4(nn_4,0) = nf5
        node_n4(nn_4,1) = n1
        node_n4(nn_4,2) = n3
        node_n4(nn_4,3) = n5
        node_n4(nn_4,4) = n7
        Grid % xn(nf5) = 0.25 * (Grid % xn(n1) + Grid % xn(n3) +  &
                                 Grid % xn(n5) + Grid % xn(n7))
        Grid % yn(nf5) = 0.25 * (Grid % yn(n1) + Grid % yn(n3) +  &
                                 Grid % yn(n5) + Grid % yn(n7))
        Grid % zn(nf5) = 0.25 * (Grid % zn(n1) + Grid % zn(n3) +  &
                                 Grid % zn(n5) + Grid % zn(n7))
      end if

      ! nf6
      do n = 1, nn_4
        if( ( node_n4(n,1) .eq. n5 .and. node_n4(n,4) .eq. n8 ) .or.  &
            ( node_n4(n,1) .eq. n8 .and. node_n4(n,4) .eq. n5 ) .or.  &
            ( node_n4(n,1) .eq. n6 .and. node_n4(n,4) .eq. n7 ) .or.  &
            ( node_n4(n,1) .eq. n7 .and. node_n4(n,4) .eq. n6 ) ) then
          nf6 = node_n4(n,0)
        end if
      end do
      if (nf6 .eq. 0) then
        nn_4 = nn_4 + 1
        Grid % n_nodes = Grid % n_nodes  + 1
        call Grid % Allocate_Nodes(Grid % n_nodes, GROWTH_MARGIN)
        nf6 = Grid % n_nodes
        node_n4(nn_4,0) = nf6
        node_n4(nn_4,1) = n5
        node_n4(nn_4,2) = n6
        node_n4(nn_4,3) = n7
        node_n4(nn_4,4) = n8
        Grid % xn(nf6) = 0.25 * (Grid % xn(n5) + Grid % xn(n6) +  &
                                 Grid % xn(n7) + Grid % xn(n8))
        Grid % yn(nf6) = 0.25 * (Grid % yn(n5) + Grid % yn(n6) +  &
                                 Grid % yn(n7) + Grid % yn(n8))
        Grid % zn(nf6) = 0.25 * (Grid % zn(n5) + Grid % zn(n6) +  &
                                 Grid % zn(n7) + Grid % zn(n8))
      end if

      !----------------------------------------!
      !   Eventually, the node in the middle   !
      !----------------------------------------!
      nn_8 = nn_8 + 1
      Grid % n_nodes = Grid % n_nodes + 1
      call Grid % Allocate_Nodes(Grid % n_nodes, GROWTH_MARGIN)
      n0  = Grid % n_nodes
      node_n8(nn_8,0) = n0
      node_n8(nn_8,1) = n1
      node_n8(nn_8,2) = n2
      node_n8(nn_8,3) = n3
      node_n8(nn_8,4) = n4
      node_n8(nn_8,5) = n5
      node_n8(nn_8,6) = n6
      node_n8(nn_8,7) = n7
      node_n8(nn_8,8) = n8
      Grid % xn(n0) = .125*(Grid % xn(n1) + Grid % xn(n2) + &
                            Grid % xn(n3) + Grid % xn(n4) + &
                            Grid % xn(n5) + Grid % xn(n6) + &
                            Grid % xn(n7) + Grid % xn(n8))
      Grid % yn(n0) = .125*(Grid % yn(n1) + Grid % yn(n2) + &
                            Grid % yn(n3) + Grid % yn(n4) + &
                            Grid % yn(n5) + Grid % yn(n6) + &
                            Grid % yn(n7) + Grid % yn(n8))
      Grid % zn(n0) = .125*(Grid % zn(n1) + Grid % zn(n2) + &
                            Grid % zn(n3) + Grid % zn(n4) + &
                            Grid % zn(n5) + Grid % zn(n6) + &
                            Grid % zn(n7) + Grid % zn(n8))

      !----------------------------!
      !   Set nodes to new cells   !
      !----------------------------!

      ! cr1 -!
      Grid % cells_n(1,cr1) = n1
      Grid % cells_n(2,cr1) = n12
      Grid % cells_n(3,cr1) = n13
      Grid % cells_n(4,cr1) = nf1
      Grid % cells_n(5,cr1) = n15
      Grid % cells_n(6,cr1) = nf2
      Grid % cells_n(7,cr1) = nf5
      Grid % cells_n(8,cr1) = n0

      ! cr2 -!
      Grid % cells_n(1,cr2) = n12
      Grid % cells_n(2,cr2) = n2
      Grid % cells_n(3,cr2) = nf1
      Grid % cells_n(4,cr2) = n24
      Grid % cells_n(5,cr2) = nf2
      Grid % cells_n(6,cr2) = n26
      Grid % cells_n(7,cr2) = n0
      Grid % cells_n(8,cr2) = nf3

      ! cr3 -!
      Grid % cells_n(1,cr3) = n13
      Grid % cells_n(2,cr3) = nf1
      Grid % cells_n(3,cr3) = n3
      Grid % cells_n(4,cr3) = n34
      Grid % cells_n(5,cr3) = nf5
      Grid % cells_n(6,cr3) = n0
      Grid % cells_n(7,cr3) = n37
      Grid % cells_n(8,cr3) = nf4

      ! cr4 -!
      Grid % cells_n(1,cr4) = nf1
      Grid % cells_n(2,cr4) = n24
      Grid % cells_n(3,cr4) = n34
      Grid % cells_n(4,cr4) = n4
      Grid % cells_n(5,cr4) = n0
      Grid % cells_n(6,cr4) = nf3
      Grid % cells_n(7,cr4) = nf4
      Grid % cells_n(8,cr4) = n48

      ! cr5 -!
      Grid % cells_n(1,cr5) = n15
      Grid % cells_n(2,cr5) = nf2
      Grid % cells_n(3,cr5) = nf5
      Grid % cells_n(4,cr5) = n0
      Grid % cells_n(5,cr5) = n5
      Grid % cells_n(6,cr5) = n56
      Grid % cells_n(7,cr5) = n57
      Grid % cells_n(8,cr5) = nf6

      ! cr6 -!
      Grid % cells_n(1,cr6) = nf2
      Grid % cells_n(2,cr6) = n26
      Grid % cells_n(3,cr6) = n0
      Grid % cells_n(4,cr6) = nf3
      Grid % cells_n(5,cr6) = n56
      Grid % cells_n(6,cr6) = n6
      Grid % cells_n(7,cr6) = nf6
      Grid % cells_n(8,cr6) = n68

      ! cr7 -!
      Grid % cells_n(1,cr7) = nf5
      Grid % cells_n(2,cr7) = n0
      Grid % cells_n(3,cr7) = n37
      Grid % cells_n(4,cr7) = nf4
      Grid % cells_n(5,cr7) = n57
      Grid % cells_n(6,cr7) = nf6
      Grid % cells_n(7,cr7) = n7
      Grid % cells_n(8,cr7) = n78

      ! cr8 -!
      Grid % cells_n(1,cr8) = n0
      Grid % cells_n(2,cr8) = nf3
      Grid % cells_n(3,cr8) = nf4
      Grid % cells_n(4,cr8) = n48
      Grid % cells_n(5,cr8) = nf6
      Grid % cells_n(6,cr8) = n68
      Grid % cells_n(7,cr8) = n78
      Grid % cells_n(8,cr8) = n8

    !-----------------------!
    !                       !
    !   Non-refined cells   !
    !                       !
    !-----------------------!
    else  ! this if started at 130

      if(ref % cell_marked(c1)) then  ! neighbor 1 refined
        Grid % cells_c( 1,c) = cell_points_to(c1)-8   &
                             + Refines_Mod_Which_Node(ref,Grid,c1,n1)
        Grid % cells_c( 7,c) = cell_points_to(c1)-8   &
                             + Refines_Mod_Which_Node(ref,Grid,c1,n2)
        Grid % cells_c(13,c) = cell_points_to(c1)-8   &
                             + Refines_Mod_Which_Node(ref,Grid,c1,n3)
        Grid % cells_c(19,c) = cell_points_to(c1)-8   &
                             + Refines_Mod_Which_Node(ref,Grid,c1,n4)
      end if

      if(ref % cell_marked(c2)) then  ! neighbor 2 refined
        Grid % cells_c( 2,c) = cell_points_to(c2)-8   &
                             + Refines_Mod_Which_Node(ref,Grid,c2,n1)
        Grid % cells_c( 8,c) = cell_points_to(c2)-8   &
                             + Refines_Mod_Which_Node(ref,Grid,c2,n2)
        Grid % cells_c(14,c) = cell_points_to(c2)-8   &
                             + Refines_Mod_Which_Node(ref,Grid,c2,n5)
        Grid % cells_c(20,c) = cell_points_to(c2)-8   &
                             + Refines_Mod_Which_Node(ref,Grid,c2,n6)
      end if

      if(ref % cell_marked(c3)) then  ! neighbor 3 refined
        Grid % cells_c( 3,c) = cell_points_to(c3)-8   &
                             + Refines_Mod_Which_Node(ref,Grid,c3,n2)
        Grid % cells_c( 9,c) = cell_points_to(c3)-8   &
                             + Refines_Mod_Which_Node(ref,Grid,c3,n4)
        Grid % cells_c(15,c) = cell_points_to(c3)-8   &
                             + Refines_Mod_Which_Node(ref,Grid,c3,n6)
        Grid % cells_c(21,c) = cell_points_to(c3)-8   &
                             + Refines_Mod_Which_Node(ref,Grid,c3,n8)
      end if

      if(ref % cell_marked(c4)) then  ! neighbor 4 refined
        Grid % cells_c( 4,c) = cell_points_to(c4)-8   &
                             + Refines_Mod_Which_Node(ref,Grid,c4,n3)
        Grid % cells_c(10,c) = cell_points_to(c4)-8   &
                             + Refines_Mod_Which_Node(ref,Grid,c4,n4)
        Grid % cells_c(16,c) = cell_points_to(c4)-8   &
                             + Refines_Mod_Which_Node(ref,Grid,c4,n7)
        Grid % cells_c(22,c) = cell_points_to(c4)-8   &
                             + Refines_Mod_Which_Node(ref,Grid,c4,n8)
      end if

      if(ref % cell_marked(c5)) then  ! neighbor 5 refined
        Grid % cells_c( 5,c) = cell_points_to(c5)-8   &
                             + Refines_Mod_Which_Node(ref,Grid,c5,n1)
        Grid % cells_c(11,c) = cell_points_to(c5)-8   &
                             + Refines_Mod_Which_Node(ref,Grid,c5,n3)
        Grid % cells_c(17,c) = cell_points_to(c5)-8   &
                             + Refines_Mod_Which_Node(ref,Grid,c5,n5)
        Grid % cells_c(23,c) = cell_points_to(c5)-8   &
                             + Refines_Mod_Which_Node(ref,Grid,c5,n7)
      end if

      if(ref % cell_marked(c6)) then  ! neighbor 6 refined
        Grid % cells_c( 6,c) = cell_points_to(c6)-8   &
                             + Refines_Mod_Which_Node(ref,Grid,c6,n5)
        Grid % cells_c(12,c) = cell_points_to(c6)-8   &
                             + Refines_Mod_Which_Node(ref,Grid,c6,n6)
        Grid % cells_c(18,c) = cell_points_to(c6)-8   &
                             + Refines_Mod_Which_Node(ref,Grid,c6,n7)
        Grid % cells_c(24,c) = cell_points_to(c6)-8   &
                             + Refines_Mod_Which_Node(ref,Grid,c6,n8)
      end if

    end if  ! if cell is refined or not
  end do    ! browse through all loops

  print *, '# Number of nodes after the refinement: ', Grid % n_nodes
  print *, '# Number of cells after the refinement: ', Grid % n_cells

  !------------------------------------------!
  !                                          !
  !   Connect the new twins, if they exist   !
  !                                          !
  !------------------------------------------!

  do na=1,nn_2
    na0=node_n2(na,0)
    na1=node_n2(na,1)
    na2=node_n2(na,2)

    if( (Grid % twin_n(na1,0) .ne. 0) .and.  &
        (Grid % twin_n(na2,0) .ne. 0) ) then

      do nb=na+1,nn_2
        nb0=node_n2(nb,0)
        nb1=node_n2(nb,1)
        nb2=node_n2(nb,2)

        if( (Grid % twin_n(nb1,0) .ne. 0) .and.  &
            (Grid % twin_n(nb2,0) .ne. 0) ) then

          if( (Grid % Are_Nodes_Twins(na1,nb1) .and.   &
               Grid % Are_Nodes_Twins(na2,nb2)) .or.   &
              (Grid % Are_Nodes_Twins(na1,nb2) .and.   &
               Grid % Are_Nodes_Twins(na2,nb1))  ) then
            if (.not. Grid % Are_Nodes_Twins(na0,nb0)) then
              Grid % twin_n(na0,0)=Grid % twin_n(na0,0)+1
              Grid % twin_n(na0,Grid % twin_n(na0,0))=nb0
              Grid % twin_n(nb0,0)=Grid % twin_n(nb0,0)+1
              Grid % twin_n(nb0,Grid % twin_n(nb0,0))=na0
            end if
          end if
        end if
      end do
    end if
  end do

  do na=1,nn_4
    na0=node_n4(na,0)
    na1=node_n4(na,1)
    na2=node_n4(na,4) ! diagonal

    if( (Grid % twin_n(na1,0) .ne. 0) .and.  &
        (Grid % twin_n(na2,0) .ne. 0) ) then

      do nb=na+1,nn_4
        nb0=node_n4(nb,0)
        nb1=node_n4(nb,1)
        nb2=node_n4(nb,4) ! diagonal

        if( (Grid % twin_n(nb1,0) .ne. 0) .and.  &
            (Grid % twin_n(nb2,0) .ne. 0) ) then

          if( (Grid % Are_Nodes_Twins(na1,nb1) .and.  &
               Grid % Are_Nodes_Twins(na2,nb2)) .or.  &
              (Grid % Are_Nodes_Twins(na1,nb2) .and.  &
               Grid % Are_Nodes_Twins(na2,nb1))  ) then
            if (.not. Grid % Are_Nodes_Twins(na0,nb0)) then
              Grid % twin_n(na0,0)=Grid % twin_n(na0,0)+1
              Grid % twin_n(na0,Grid % twin_n(na0,0))=nb0
              Grid % twin_n(nb0,0)=Grid % twin_n(nb0,0)+1
              Grid % twin_n(nb0,Grid % twin_n(nb0,0))=na0
            end if
          end if
        end if
      end do
    end if
  end do

  !----------------------------!
  !                            !
  !   Delete redundant cells   !
  !                            !
  !----------------------------!

  ! Initialize the new numbers for the cells
  do c = -Grid % n_bnd_cells, Grid % n_cells
    Grid % new_c(c) = c
  end do

  del = 0
  do c = 1, Grid % n_cells
    if(ref % cell_marked(c)) then
      Grid % new_c(c) = -1
      del = del+1
    else
      Grid % new_c(c) = c - del
    end if
  end do
  print *, '# Deleted cells:', del

  do c = 1, Grid % n_cells
    if(Grid % new_c(c) .ne. -1) then

      ! Update the cell numbers. Watch out ! The numbers you are
      ! updating are old, so double indexing is needed
      do n = 1, 24  ! n is neighbour now
        Grid % cells_c( n, Grid % new_c(c) )  &
          = Grid % new_c( Grid % cells_c( n, c ) )
      end do

      ! Update the node numbers
      do n = 1, 8   ! n is node now
        Grid % cells_n( n, Grid % new_c(c) ) = Grid % cells_n( n, c )
      end do

      ! The line below was never checked !
      ref % cell_level(Grid % new_c(c)) = ref % cell_level(c)
    end if
  end do

  do c=Grid % n_cells-del+1, Grid % n_cells       ! erase old data
    do n = 1, 24                                  ! n is neighbour now
      Grid % cells_c( n, c ) = 0
    end do
  end do

  Grid % n_cells = Grid % n_cells - del

  print *, '# Number of cells after the renumeration: ', Grid % n_cells

  end subroutine
