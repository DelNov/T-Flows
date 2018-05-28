!==============================================================================!
  subroutine Refine_Marked_Cells(grid, lev)
!------------------------------------------------------------------------------!
!   Refine the marked cells.                                                   !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use gen_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: lev 
!----------------------------------[Calling]-----------------------------------!
  logical :: Are_Nodes_Twins
  integer :: Which_Node
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, n_cells_old, c1, c2, c3, c4, c5, c6
  integer :: cr1, cr2, cr3, cr4, cr5, cr6, cr7, cr8
  integer :: n, n_nodes_old, n1, n2, n3, n4, n5, n6, n7, n8
  integer :: n12,n13,n24,n34,n15,n26,n37,n48,n56,n57,n68,n78
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
  print *, '# Refine: Number of nodes: ', grid % n_nodes 
  print *, '#         Number of cells: ', grid % n_cells 
  print *, '#----------------------------------------------'

  n_cells_old = grid % n_cells 
  n_nodes_old = grid % n_nodes
  nn_2 = 0 
  nn_4 = 0
  nn_8 = 0

  ! Allocate local memory
  allocate (node_n2(grid % max_n_nodes,0:2));     node_n2        = 0 
  allocate (node_n4(grid % max_n_nodes,0:4));     node_n4        = 0
  allocate (node_n8(grid % max_n_nodes,0:8));     node_n8        = 0
  allocate (cell_points_to(grid % max_n_nodes));  cell_points_to = 0

  !---------------------!
  !                     !
  !   Count new celss   !
  !                     !
  !---------------------!
  do c = 1, n_cells_old
    if(cell_marked(c)) then 
      grid % n_cells = grid % n_cells + 8
      cell_points_to(c) = grid % n_cells   ! now points to cr8
    end if
  end do

  do c = 1, n_cells_old

    c1 = grid % cells_c(1,c)
    c2 = grid % cells_c(2,c)
    c3 = grid % cells_c(3,c)
    c4 = grid % cells_c(4,c)
    c5 = grid % cells_c(5,c)
    c6 = grid % cells_c(6,c)

    n1 = grid % cells_n(1,c)
    n2 = grid % cells_n(2,c)
    n3 = grid % cells_n(3,c)
    n4 = grid % cells_n(4,c)
    n5 = grid % cells_n(5,c)
    n6 = grid % cells_n(6,c)
    n7 = grid % cells_n(7,c)
    n8 = grid % cells_n(8,c)

    !-------------------!
    !                   !
    !   Refined cells   !
    !                   !
    !-------------------!
    if( cell_marked(c) ) then  ! only refined

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

      grid % material(cr1) = grid % material(c) 
      grid % material(cr2) = grid % material(c) 
      grid % material(cr3) = grid % material(c) 
      grid % material(cr4) = grid % material(c) 
      grid % material(cr5) = grid % material(c) 
      grid % material(cr6) = grid % material(c) 
      grid % material(cr7) = grid % material(c) 
      grid % material(cr8) = grid % material(c) 

      !-----------------------------------------------!
      !   Internal links do not depend on neighbors   !
      !-----------------------------------------------!

      ! 6
      grid % cells_c(6,cr1) = cr5
      grid % cells_c(6,cr2) = cr6
      grid % cells_c(6,cr3) = cr7
      grid % cells_c(6,cr4) = cr8

      ! 5
      grid % cells_c(5,cr2) = cr1
      grid % cells_c(5,cr4) = cr3
      grid % cells_c(5,cr6) = cr5
      grid % cells_c(5,cr8) = cr7

      ! 4
      grid % cells_c(4,cr1) = cr3
      grid % cells_c(4,cr2) = cr4
      grid % cells_c(4,cr5) = cr7
      grid % cells_c(4,cr6) = cr8

      ! 3
      grid % cells_c(3,cr1) = cr2
      grid % cells_c(3,cr3) = cr4
      grid % cells_c(3,cr5) = cr6
      grid % cells_c(3,cr7) = cr8

      ! 2
      grid % cells_c(2,cr3) = cr1
      grid % cells_c(2,cr4) = cr2
      grid % cells_c(2,cr7) = cr5
      grid % cells_c(2,cr8) = cr6

      ! 1
      grid % cells_c(1,cr5) = cr1
      grid % cells_c(1,cr6) = cr2
      grid % cells_c(1,cr7) = cr3
      grid % cells_c(1,cr8) = cr4

      !-------------------------!
      !   Level of refinement   !
      !-------------------------!
      level(cr1) = lev
      level(cr2) = lev
      level(cr3) = lev
      level(cr4) = lev
      level(cr5) = lev
      level(cr6) = lev
      level(cr7) = lev
      level(cr8) = lev

      !----------------------------------------!         
      !   External links depend on neighbors   !
      !----------------------------------------!
      if(.not. cell_marked(c1)) then  ! neighbor 1 not refined
        grid % cells_c(1,cr1) = c1
        grid % cells_c(1,cr2) = c1
        grid % cells_c(1,cr3) = c1
        grid % cells_c(1,cr4) = c1
      else                        ! neighbor 1 refined
        grid % cells_c(1,cr1) = cell_points_to(c1) - 8 + Which_Node(grid, c1,n1)
        grid % cells_c(1,cr2) = cell_points_to(c1) - 8 + Which_Node(grid, c1,n2)
        grid % cells_c(1,cr3) = cell_points_to(c1) - 8 + Which_Node(grid, c1,n3)
        grid % cells_c(1,cr4) = cell_points_to(c1) - 8 + Which_Node(grid, c1,n4)
      endif           

      if(.not. cell_marked(c2)) then  ! neighbor 2 not refined
        grid % cells_c(2,cr1) = c2
        grid % cells_c(2,cr2) = c2
        grid % cells_c(2,cr5) = c2
        grid % cells_c(2,cr6) = c2
      else                        ! neighbor 2 refined
        grid % cells_c(2,cr1) = cell_points_to(c2) - 8 + Which_Node(grid, c2, n1)
        grid % cells_c(2,cr2) = cell_points_to(c2) - 8 + Which_Node(grid, c2, n2)
        grid % cells_c(2,cr5) = cell_points_to(c2) - 8 + Which_Node(grid, c2, n5)
        grid % cells_c(2,cr6) = cell_points_to(c2) - 8 + Which_Node(grid, c2, n6)
      endif           

      if(.not. cell_marked(c3)) then  ! neighbor 3 not refined
        grid % cells_c(3,cr2) = c3
        grid % cells_c(3,cr4) = c3
        grid % cells_c(3,cr6) = c3
        grid % cells_c(3,cr8) = c3
      else                        ! neighbor 3 refined
        grid % cells_c(3,cr2) = cell_points_to(c3) - 8 + Which_Node(grid, c3, n2)
        grid % cells_c(3,cr4) = cell_points_to(c3) - 8 + Which_Node(grid, c3, n4)
        grid % cells_c(3,cr6) = cell_points_to(c3) - 8 + Which_Node(grid, c3, n6)
        grid % cells_c(3,cr8) = cell_points_to(c3) - 8 + Which_Node(grid, c3, n8)
      endif           

      if(.not. cell_marked(c4)) then  ! neighbor 4 not refine
        grid % cells_c(4,cr3) = c4
        grid % cells_c(4,cr4) = c4
        grid % cells_c(4,cr7) = c4
        grid % cells_c(4,cr8) = c4
      else                        ! neighbor 4 refine
        grid % cells_c(4,cr3) = cell_points_to(c4) - 8 + Which_Node(grid, c4, n3)
        grid % cells_c(4,cr4) = cell_points_to(c4) - 8 + Which_Node(grid, c4, n4)
        grid % cells_c(4,cr7) = cell_points_to(c4) - 8 + Which_Node(grid, c4, n7)
        grid % cells_c(4,cr8) = cell_points_to(c4) - 8 + Which_Node(grid, c4, n8)
      endif           

      if(.not. cell_marked(c5)) then  ! neighbor 5 not refined
        grid % cells_c(5,cr1) = c5
        grid % cells_c(5,cr3) = c5
        grid % cells_c(5,cr5) = c5
        grid % cells_c(5,cr7) = c5
      else                        ! neighbor 5 refined
        grid % cells_c(5,cr1) = cell_points_to(c5) - 8 + Which_Node(grid, c5, n1)
        grid % cells_c(5,cr3) = cell_points_to(c5) - 8 + Which_Node(grid, c5, n3)
        grid % cells_c(5,cr5) = cell_points_to(c5) - 8 + Which_Node(grid, c5, n5)
        grid % cells_c(5,cr7) = cell_points_to(c5) - 8 + Which_Node(grid, c5, n7)
      endif

      if(.not. cell_marked(c6)) then  ! neighbor 6 not refined
        grid % cells_c(6,cr5) = c6
        grid % cells_c(6,cr6) = c6
        grid % cells_c(6,cr7) = c6
        grid % cells_c(6,cr8) = c6
      else                        ! neighbor 6 refined
        grid % cells_c(6,cr5) = cell_points_to(c6) - 8 + Which_Node(grid, c6, n5)
        grid % cells_c(6,cr6) = cell_points_to(c6) - 8 + Which_Node(grid, c6, n6)
        grid % cells_c(6,cr7) = cell_points_to(c6) - 8 + Which_Node(grid, c6, n7)
        grid % cells_c(6,cr8) = cell_points_to(c6) - 8 + Which_Node(grid, c6, n8)
      endif           

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
      grid % n_nodes  = grid % n_nodes  + 1
      n12 = grid % n_nodes
      node_n2(nn_2,0) = n12
      node_n2(nn_2,1) = n1
      node_n2(nn_2,2) = n2
      grid % xn(n12) = .5 * (grid % xn(n1) + grid % xn(n2))
      grid % yn(n12) = .5 * (grid % yn(n1) + grid % yn(n2))
      grid % zn(n12) = .5 * (grid % zn(n1) + grid % zn(n2))
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
      grid % n_nodes  = grid % n_nodes  + 1
      n13 = grid % n_nodes
      node_n2(nn_2,0) = n13
      node_n2(nn_2,1) = n1
      node_n2(nn_2,2) = n3
      grid % xn(n13) = .5 * (grid % xn(n1) + grid % xn(n3))
      grid % yn(n13) = .5 * (grid % yn(n1) + grid % yn(n3))
      grid % zn(n13) = .5 * (grid % zn(n1) + grid % zn(n3))
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
      grid % n_nodes  = grid % n_nodes  + 1
      n24 = grid % n_nodes
      node_n2(nn_2,0) = n24
      node_n2(nn_2,1) = n2
      node_n2(nn_2,2) = n4
      grid % xn(n24) = .5 * (grid % xn(n2) + grid % xn(n4))
      grid % yn(n24) = .5 * (grid % yn(n2) + grid % yn(n4))
      grid % zn(n24) = .5 * (grid % zn(n2) + grid % zn(n4))
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
      grid % n_nodes  = grid % n_nodes  + 1
      n34 = grid % n_nodes
      node_n2(nn_2,0) = n34
      node_n2(nn_2,1) = n3
      node_n2(nn_2,2) = n4
      grid % xn(n34) = .5 * (grid % xn(n3) + grid % xn(n4))
      grid % yn(n34) = .5 * (grid % yn(n3) + grid % yn(n4))
      grid % zn(n34) = .5 * (grid % zn(n3) + grid % zn(n4))
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
      grid % n_nodes  = grid % n_nodes  + 1
      n15 = grid % n_nodes
      node_n2(nn_2,0) = n15
      node_n2(nn_2,1) = n1
      node_n2(nn_2,2) = n5
      grid % xn(n15) = .5 * (grid % xn(n1) + grid % xn(n5))
      grid % yn(n15) = .5 * (grid % yn(n1) + grid % yn(n5))
      grid % zn(n15) = .5 * (grid % zn(n1) + grid % zn(n5))
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
      grid % n_nodes  = grid % n_nodes  + 1
      n26 = grid % n_nodes
      node_n2(nn_2,0) = n26
      node_n2(nn_2,1) = n2
      node_n2(nn_2,2) = n6
      grid % xn(n26) = .5 * (grid % xn(n2) + grid % xn(n6))
      grid % yn(n26) = .5 * (grid % yn(n2) + grid % yn(n6))
      grid % zn(n26) = .5 * (grid % zn(n2) + grid % zn(n6))
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
      grid % n_nodes  = grid % n_nodes  + 1
      n37 = grid % n_nodes
      node_n2(nn_2,0) = n37
      node_n2(nn_2,1) = n3
      node_n2(nn_2,2) = n7
      grid % xn(n37) = .5 * (grid % xn(n3) + grid % xn(n7))
      grid % yn(n37) = .5 * (grid % yn(n3) + grid % yn(n7))
      grid % zn(n37) = .5 * (grid % zn(n3) + grid % zn(n7))
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
      grid % n_nodes  = grid % n_nodes  + 1
      n48 = grid % n_nodes
      node_n2(nn_2,0) = n48
      node_n2(nn_2,1) = n4
      node_n2(nn_2,2) = n8
      grid % xn(n48) = .5 * (grid % xn(n4) + grid % xn(n8))
      grid % yn(n48) = .5 * (grid % yn(n4) + grid % yn(n8))
      grid % zn(n48) = .5 * (grid % zn(n4) + grid % zn(n8))
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
      grid % n_nodes  = grid % n_nodes  + 1
      n56 = grid % n_nodes
      node_n2(nn_2,0) = n56
      node_n2(nn_2,1) = n5
      node_n2(nn_2,2) = n6
      grid % xn(n56) = .5 * (grid % xn(n5) + grid % xn(n6))
      grid % yn(n56) = .5 * (grid % yn(n5) + grid % yn(n6))
      grid % zn(n56) = .5 * (grid % zn(n5) + grid % zn(n6))
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
      grid % n_nodes  = grid % n_nodes  + 1
      n57 = grid % n_nodes
      node_n2(nn_2,0) = n57
      node_n2(nn_2,1) = n5
      node_n2(nn_2,2) = n7
      grid % xn(n57) = .5 * (grid % xn(n5) + grid % xn(n7))
      grid % yn(n57) = .5 * (grid % yn(n5) + grid % yn(n7))
      grid % zn(n57) = .5 * (grid % zn(n5) + grid % zn(n7))
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
      grid % n_nodes  = grid % n_nodes  + 1
      n68 = grid % n_nodes
      node_n2(nn_2,0) = n68
      node_n2(nn_2,1) = n6
      node_n2(nn_2,2) = n8
      grid % xn(n68) = .5 * (grid % xn(n6) + grid % xn(n8))
      grid % yn(n68) = .5 * (grid % yn(n6) + grid % yn(n8))
      grid % zn(n68) = .5 * (grid % zn(n6) + grid % zn(n8))
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
      grid % n_nodes  = grid % n_nodes  + 1
      n78 = grid % n_nodes
      node_n2(nn_2,0) = n78
      node_n2(nn_2,1) = n7
      node_n2(nn_2,2) = n8
      grid % xn(n78) = .5 * (grid % xn(n7) + grid % xn(n8))
      grid % yn(n78) = .5 * (grid % yn(n7) + grid % yn(n8))
      grid % zn(n78) = .5 * (grid % zn(n7) + grid % zn(n8))
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
      grid % n_nodes  = grid % n_nodes  + 1
      nf1 = grid % n_nodes
      node_n4(nn_4,0) = nf1
      node_n4(nn_4,1) = n1
      node_n4(nn_4,2) = n2
      node_n4(nn_4,3) = n3
      node_n4(nn_4,4) = n4
      grid % xn(nf1) = 0.25 * (grid % xn(n1) + grid % xn(n2) +  &
                               grid % xn(n3) + grid % xn(n4))
      grid % yn(nf1) = 0.25 * (grid % yn(n1) + grid % yn(n2) +  &
                               grid % yn(n3) + grid % yn(n4))
      grid % zn(nf1) = 0.25 * (grid % zn(n1) + grid % zn(n2) +  &
                               grid % zn(n3) + grid % zn(n4))
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
      grid % n_nodes  = grid % n_nodes  + 1
      nf2 = grid % n_nodes
      node_n4(nn_4,0) = nf2
      node_n4(nn_4,1) = n1
      node_n4(nn_4,2) = n2
      node_n4(nn_4,3) = n5
      node_n4(nn_4,4) = n6
      grid % xn(nf2) = 0.25 * (grid % xn(n1) + grid % xn(n2) +  &
                               grid % xn(n5) + grid % xn(n6))
      grid % yn(nf2) = 0.25 * (grid % yn(n1) + grid % yn(n2) +  &
                               grid % yn(n5) + grid % yn(n6))
      grid % zn(nf2) = 0.25 * (grid % zn(n1) + grid % zn(n2) +  &
                               grid % zn(n5) + grid % zn(n6))
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
      grid % n_nodes  = grid % n_nodes  + 1
      nf3 = grid % n_nodes
      node_n4(nn_4,0) = nf3
      node_n4(nn_4,1) = n2
      node_n4(nn_4,2) = n4
      node_n4(nn_4,3) = n6
      node_n4(nn_4,4) = n8
      grid % xn(nf3) = 0.25 * (grid % xn(n2) + grid % xn(n4) +  &
                               grid % xn(n6) + grid % xn(n8))
      grid % yn(nf3) = 0.25 * (grid % yn(n2) + grid % yn(n4) +  &
                               grid % yn(n6) + grid % yn(n8))
      grid % zn(nf3) = 0.25 * (grid % zn(n2) + grid % zn(n4) +  &
                               grid % zn(n6) + grid % zn(n8))
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
      grid % n_nodes  = grid % n_nodes  + 1
      nf4 = grid % n_nodes
      node_n4(nn_4,0) = nf4
      node_n4(nn_4,1) = n3
      node_n4(nn_4,2) = n4
      node_n4(nn_4,3) = n7
      node_n4(nn_4,4) = n8
      grid % xn(nf4) = 0.25 * (grid % xn(n3) + grid % xn(n4) +  &
                               grid % xn(n7) + grid % xn(n8))
      grid % yn(nf4) = 0.25 * (grid % yn(n3) + grid % yn(n4) +  &
                               grid % yn(n7) + grid % yn(n8))
      grid % zn(nf4) = 0.25 * (grid % zn(n3) + grid % zn(n4) +  &
                               grid % zn(n7) + grid % zn(n8))
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
      grid % n_nodes  = grid % n_nodes  + 1
      nf5 = grid % n_nodes
      node_n4(nn_4,0) = nf5
      node_n4(nn_4,1) = n1
      node_n4(nn_4,2) = n3
      node_n4(nn_4,3) = n5
      node_n4(nn_4,4) = n7
      grid % xn(nf5) = 0.25 * (grid % xn(n1) + grid % xn(n3) +  &
                               grid % xn(n5) + grid % xn(n7))
      grid % yn(nf5) = 0.25 * (grid % yn(n1) + grid % yn(n3) +  &
                               grid % yn(n5) + grid % yn(n7))
      grid % zn(nf5) = 0.25 * (grid % zn(n1) + grid % zn(n3) +  &
                               grid % zn(n5) + grid % zn(n7))
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
      grid % n_nodes  = grid % n_nodes  + 1
      nf6 = grid % n_nodes
      node_n4(nn_4,0) = nf6
      node_n4(nn_4,1) = n5
      node_n4(nn_4,2) = n6
      node_n4(nn_4,3) = n7
      node_n4(nn_4,4) = n8
      grid % xn(nf6) = 0.25 * (grid % xn(n5) + grid % xn(n6) +  &
                               grid % xn(n7) + grid % xn(n8))
      grid % yn(nf6) = 0.25 * (grid % yn(n5) + grid % yn(n6) +  &
                               grid % yn(n7) + grid % yn(n8))
      grid % zn(nf6) = 0.25 * (grid % zn(n5) + grid % zn(n6) +  &
                               grid % zn(n7) + grid % zn(n8))
    end if 

    !----------------------------------------!
    !   Eventually, the node in the middle   !
    !----------------------------------------!
    nn_8 = nn_8 + 1
    grid % n_nodes  = grid % n_nodes + 1
    n0  = grid % n_nodes
    node_n8(nn_8,0) = n0 
    node_n8(nn_8,1) = n1
    node_n8(nn_8,2) = n2
    node_n8(nn_8,3) = n3
    node_n8(nn_8,4) = n4
    node_n8(nn_8,5) = n5
    node_n8(nn_8,6) = n6
    node_n8(nn_8,7) = n7
    node_n8(nn_8,8) = n8
    grid % xn(n0) = .125*(grid % xn(n1) + grid % xn(n2) + &
                          grid % xn(n3) + grid % xn(n4) + &
                          grid % xn(n5) + grid % xn(n6) + &
                          grid % xn(n7) + grid % xn(n8))
    grid % yn(n0) = .125*(grid % yn(n1) + grid % yn(n2) + &
                          grid % yn(n3) + grid % yn(n4) + &
                          grid % yn(n5) + grid % yn(n6) + &
                          grid % yn(n7) + grid % yn(n8))
    grid % zn(n0) = .125*(grid % zn(n1) + grid % zn(n2) + &
                          grid % zn(n3) + grid % zn(n4) + &
                          grid % zn(n5) + grid % zn(n6) + &
                          grid % zn(n7) + grid % zn(n8))

    !----------------------------!
    !   Set nodes to new cells   !
    !----------------------------!

    ! cr1 -!
    grid % cells_n(1,cr1) = n1
    grid % cells_n(2,cr1) = n12
    grid % cells_n(3,cr1) = n13
    grid % cells_n(4,cr1) = nf1
    grid % cells_n(5,cr1) = n15
    grid % cells_n(6,cr1) = nf2
    grid % cells_n(7,cr1) = nf5
    grid % cells_n(8,cr1) = n0 

    ! cr2 -!
    grid % cells_n(1,cr2) = n12
    grid % cells_n(2,cr2) = n2 
    grid % cells_n(3,cr2) = nf1
    grid % cells_n(4,cr2) = n24
    grid % cells_n(5,cr2) = nf2
    grid % cells_n(6,cr2) = n26 
    grid % cells_n(7,cr2) = n0 
    grid % cells_n(8,cr2) = nf3

    ! cr3 -!
    grid % cells_n(1,cr3) = n13
    grid % cells_n(2,cr3) = nf1
    grid % cells_n(3,cr3) = n3 
    grid % cells_n(4,cr3) = n34
    grid % cells_n(5,cr3) = nf5
    grid % cells_n(6,cr3) = n0  
    grid % cells_n(7,cr3) = n37
    grid % cells_n(8,cr3) = nf4

    ! cr4 -!
    grid % cells_n(1,cr4) = nf1
    grid % cells_n(2,cr4) = n24
    grid % cells_n(3,cr4) = n34
    grid % cells_n(4,cr4) = n4 
    grid % cells_n(5,cr4) = n0 
    grid % cells_n(6,cr4) = nf3 
    grid % cells_n(7,cr4) = nf4
    grid % cells_n(8,cr4) = n48

    ! cr5 -!
    grid % cells_n(1,cr5) = n15
    grid % cells_n(2,cr5) = nf2
    grid % cells_n(3,cr5) = nf5
    grid % cells_n(4,cr5) = n0 
    grid % cells_n(5,cr5) = n5 
    grid % cells_n(6,cr5) = n56 
    grid % cells_n(7,cr5) = n57
    grid % cells_n(8,cr5) = nf6

    ! cr6 -!
    grid % cells_n(1,cr6) = nf2
    grid % cells_n(2,cr6) = n26
    grid % cells_n(3,cr6) = n0 
    grid % cells_n(4,cr6) = nf3
    grid % cells_n(5,cr6) = n56
    grid % cells_n(6,cr6) = n6  
    grid % cells_n(7,cr6) = nf6
    grid % cells_n(8,cr6) = n68

    ! cr7 -!
    grid % cells_n(1,cr7) = nf5
    grid % cells_n(2,cr7) = n0 
    grid % cells_n(3,cr7) = n37
    grid % cells_n(4,cr7) = nf4
    grid % cells_n(5,cr7) = n57
    grid % cells_n(6,cr7) = nf6 
    grid % cells_n(7,cr7) = n7 
    grid % cells_n(8,cr7) = n78

    ! cr8 -!
    grid % cells_n(1,cr8) = n0 
    grid % cells_n(2,cr8) = nf3
    grid % cells_n(3,cr8) = nf4
    grid % cells_n(4,cr8) = n48
    grid % cells_n(5,cr8) = nf6
    grid % cells_n(6,cr8) = n68 
    grid % cells_n(7,cr8) = n78
    grid % cells_n(8,cr8) = n8 

    !-----------------------!
    !                       !
    !   Non-refined cells   !
    !                       !
    !-----------------------!
    else

      if(cell_marked(c1)) then  ! neighbor 1 refined
        grid % cells_c( 1,c) = cell_points_to(c1) - 8 + Which_Node(grid, c1,n1)
        grid % cells_c( 7,c) = cell_points_to(c1) - 8 + Which_Node(grid, c1,n2)
        grid % cells_c(13,c) = cell_points_to(c1) - 8 + Which_Node(grid, c1,n3)
        grid % cells_c(19,c) = cell_points_to(c1) - 8 + Which_Node(grid, c1,n4)
      endif         

      if(cell_marked(c2)) then  ! neighbor 2 refined
        grid % cells_c( 2,c) = cell_points_to(c2) - 8 + Which_Node(grid, c2, n1)
        grid % cells_c( 8,c) = cell_points_to(c2) - 8 + Which_Node(grid, c2, n2)
        grid % cells_c(14,c) = cell_points_to(c2) - 8 + Which_Node(grid, c2, n5)
        grid % cells_c(20,c) = cell_points_to(c2) - 8 + Which_Node(grid, c2, n6)
      endif           

      if(cell_marked(c3)) then  ! neighbor 3 refined
        grid % cells_c( 3,c) = cell_points_to(c3) - 8 + Which_Node(grid, c3, n2)
        grid % cells_c( 9,c) = cell_points_to(c3) - 8 + Which_Node(grid, c3, n4)
        grid % cells_c(15,c) = cell_points_to(c3) - 8 + Which_Node(grid, c3, n6)
        grid % cells_c(21,c) = cell_points_to(c3) - 8 + Which_Node(grid, c3, n8)
      endif           

      if(cell_marked(c4)) then  ! neighbor 4 refined
        grid % cells_c( 4,c) = cell_points_to(c4) - 8 + Which_Node(grid, c4, n3)
        grid % cells_c(10,c) = cell_points_to(c4) - 8 + Which_Node(grid, c4, n4)
        grid % cells_c(16,c) = cell_points_to(c4) - 8 + Which_Node(grid, c4, n7)
        grid % cells_c(22,c) = cell_points_to(c4) - 8 + Which_Node(grid, c4, n8)
      endif           

      if(cell_marked(c5)) then  ! neighbor 5 refined
        grid % cells_c( 5,c) = cell_points_to(c5) - 8 + Which_Node(grid, c5, n1)
        grid % cells_c(11,c) = cell_points_to(c5) - 8 + Which_Node(grid, c5, n3)
        grid % cells_c(17,c) = cell_points_to(c5) - 8 + Which_Node(grid, c5, n5)
        grid % cells_c(23,c) = cell_points_to(c5) - 8 + Which_Node(grid, c5, n7)
      endif

      if(cell_marked(c6)) then  ! neighbor 6 refined
        grid % cells_c( 6,c) = cell_points_to(c6) - 8 + Which_Node(grid, c6, n5)
        grid % cells_c(12,c) = cell_points_to(c6) - 8 + Which_Node(grid, c6, n6)
        grid % cells_c(18,c) = cell_points_to(c6) - 8 + Which_Node(grid, c6, n7)
        grid % cells_c(24,c) = cell_points_to(c6) - 8 + Which_Node(grid, c6, n8)
      endif           

    end if   
  end do

  print *, '# Number of nodes after the refinement: ', grid % n_nodes 
  print *, '# Number of cells after the refinement: ', grid % n_cells 

  !------------------------------------------!
  !                                          !
  !   Connect the new twins, if they exist   !
  !                                          !
  !------------------------------------------!  

  do na=1,nn_2
    na0=node_n2(na,0)
    na1=node_n2(na,1)
    na2=node_n2(na,2)

    if( (twin_n(na1,0) .ne. 0).and.(twin_n(na2,0) .ne. 0) ) then

      do nb=na+1,nn_2
        nb0=node_n2(nb,0)
        nb1=node_n2(nb,1)
        nb2=node_n2(nb,2)

        if( (twin_n(nb1,0) .ne. 0).and.(twin_n(nb2,0) .ne. 0) ) then

          if( (Are_Nodes_Twins(na1,nb1) .and.   &
               Are_Nodes_Twins(na2,nb2)) .or.   &
              (Are_Nodes_Twins(na1,nb2) .and.   &
               Are_Nodes_Twins(na2,nb1))  ) then
            if (.not. Are_Nodes_Twins(na0,nb0)) then
              twin_n(na0,0)=twin_n(na0,0)+1
              twin_n(na0,twin_n(na0,0))=nb0
              twin_n(nb0,0)=twin_n(nb0,0)+1
              twin_n(nb0,twin_n(nb0,0))=na0
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

    if( (twin_n(na1,0) .ne. 0).and.(twin_n(na2,0) .ne. 0) ) then

      do nb=na+1,nn_4
        nb0=node_n4(nb,0)
        nb1=node_n4(nb,1)
        nb2=node_n4(nb,4) ! diagonal

        if( (twin_n(nb1,0) .ne. 0).and.(twin_n(nb2,0) .ne. 0) ) then

          if( (Are_Nodes_Twins(na1,nb1) .and. Are_Nodes_Twins(na2,nb2)) .or.          &
              (Are_Nodes_Twins(na1,nb2) .and. Are_Nodes_Twins(na2,nb1))  ) then
            if (.not. Are_Nodes_Twins(na0,nb0)) then
              twin_n(na0,0)=twin_n(na0,0)+1
              twin_n(na0,twin_n(na0,0))=nb0
              twin_n(nb0,0)=twin_n(nb0,0)+1
              twin_n(nb0,twin_n(nb0,0))=na0
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
  do c = -grid % n_bnd_cells, grid % n_cells
    new_n(c)=c
  end do

  del = 0 
  do c = 1, grid % n_cells
    if(cell_marked(c)) then
      new_n(c) = -1
      del = del+1
    else
      new_n(c) = c - del 
    endif 
  end do
  print *, '# Deleted cells:', del

  do c = 1, grid % n_cells
    if(new_n(c) .ne. -1) then

      ! Update the cell numbers. Watch out ! The numbers you are
      ! updating are old, so double indexing is needed
      do n = 1, 24  ! n is neighbour now
        grid % cells_c( n, new_n(c) ) = new_n( grid % cells_c( n, c ) )
      end do

      ! Update the node numbers
      do n = 1, 8   ! n is node now
        grid % cells_n( n, new_n(c) ) = grid % cells_n( n, c )
      end do

      grid % material( new_n(c) ) = grid % material( c )  ! -> never checked !
      level( new_n(c) )    = level( c )                   ! -> never checked !
    end if
  end do

  do c=grid % n_cells-del+1, grid % max_n_nodes   ! erase old data
    do n = 1, 24                                  ! n is neighbour now
      grid % cells_c( n, c ) = 0
    end do
  end do

  grid % n_cells = grid % n_cells - del    

  print *, '# Number of cells after the renumeration: ', grid % n_cells 

  deallocate(node_n2)
  deallocate(node_n4)
  deallocate(node_n8)
  deallocate(cell_points_to)

  end subroutine
