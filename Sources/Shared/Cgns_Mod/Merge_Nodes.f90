!==============================================================================!
  subroutine Cgns_Mod_Merge_Nodes(grid)
!------------------------------------------------------------------------------!
!   For each interface in geometry merges nodes on interfaces and remaps       !
!   cell_connections                                                           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
  use Sort_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!------------------------------------------------------------------------------!
  include '../Shared/Approx_Real.int'
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, cnt_node, n, i, k, v
  real,    allocatable :: criterion(:) ! sorting criterion
  integer, allocatable :: old_seq(:), new_seq(:)
  real,    allocatable :: x_new(:), y_new(:), z_new(:)
  real                 :: big, small
  integer              :: int
  integer              :: cnt_nodes_on_int_to_keep, cnt_nodes_on_int_to_rem
  integer              :: cnt_nodes_on_int_total
  integer              :: n1, n2
  logical, allocatable :: nodes_to_remove(:) ! marked duplicated nodes to remove

!==============================================================================!

  print '(a)', ' #=============================================='
  print '(a)', ' # Merging blocks since they have common interfaces'
  print '(a)', ' # Hint: Join blocks in mesher to avoid issues'
  print '(a,i26)', ' # Old number of nodes: ', grid % n_nodes
  print '(a)', ' #----------------------------------------------'
  ! '/2' because mixed interface type had weight 1, while quad and tri had 2
  cnt_int = cnt_int / 2

  !----------------------------------------------------------------------------!
  !   At this point number of unique interfaces cnt_int is known               !
  !   Cells on interfaces are in interface_cells                               !
  !   Array interface_cells has 4 indices: (1:2, cnt_int_cells, 1:4, cnt_int)  !
  !                                                                            !
  !   Each interface has two domains from 2 different blocks                   !
  !   First index denotes this first and second part of a interface            !
  !   Second index is interface cell id                                        !
  !   Third index is node id for a cell. -1 mean this node is not specified    !
  !   Last index tells to which unique interface this cell belongs.            !
  !                                                                            !
  !   Unfortunately coordinates of cell nodes inside domains                   !
  !   interface_cells(1,...) and interface_cells(2,...)                        !
  !   (interface_nodes(n, int) = 1 and interface_nodes(n, int) = -1)           !
  !   do not match, because they originated from different blocks.             !
  !   That is why sorting is required by some criterion to group them.         !
  !----------------------------------------------------------------------------!
  !   Original block structure with duplicate nodes:                           !
  !        x  y  z                                                             !
  !    1:  1  1  1         13--------9    <->      5--------1                  !
  !    2:  1  1  0        /|        /|    <->     /|       /|                  !
  !    3:  1 -1  1       / |       / |    <->    / |      / |                  !
  !    4:  1 -1  0      14-|------10 |    <->   6--|-----2--|                  !
  !    5:  0  1  1      |  15     |  11   <->   |  7     |  3                  !
  !    6:  0  1  0      |  /      | /     <->   | /      | /                   !
  !    7:  0 -1  1      | /       |/      <->   |/       |/                    !
  !    8:  0 -1  0      16--------12      <->   8--------4                     !
  !    9:  0  1  1                                                             !
  !   10:  0  1  0        block 2                    block 1                   !
  !   11:  0 -1  1    y                                                        !
  !   12:  0 -1  0    ^                                                        !
  !   13: -1  1  1    |   ^ z                                                  !
  !   14: -1  1  0    |  /                                                     !
  !   15: -1 -1  1    | /                                                      !
  !   16: -1 -1  0    --------> x                                              !
  !                                                                            !
  !   block 1: 1   2   3   4   5   6   7   8                                   !
  !   block 2: 9  10  11  12  13  14  15  16                                   !
  !                                                                            !
  !   after function:                                                          !
  !        x  y  z                                                             !
  !    1:  1  1  1       9----------5   <->      5--------1                    !
  !    2:  1  1  0      / |        /|   <->     /|       /|                    !
  !    3:  1 -1  1     /  |       / |   <->    / |      / |                    !
  !    4:  1 -1  0    10--|------6--|   <->   6--|-----2--|                    !
  !    5:  0  1  1    |   11     |  7   <->   |  7     |  3                    !
  !    6:  0  1  0    |  /       | /    <->   | /      | /                     !
  !    7:  0 -1  1    | /        |/     <->   |/       |/                      !
  !    8:  0 -1  0    12---------8      <->   8--------4                       !
  !    9: -1  1  1                                                             !
  !   10: -1  1  0    block 2                    block 1                       !
  !   11: -1 -1  1                                                             !
  !   12: -1 -1  0                                                             !
  !                                                                            !
  !   block 1: 1  2  3  4   5   6   7   8                                      !
  !   block 2: 5  6  7  8   9  10  11  12                                      !
  !                                                                            !
  !----------------------------------------------------------------------------!

  ! Estimate big and small
  call Grid_Mod_Estimate_Big_And_Small(grid, big, small)

  if (verbose) then
    print '(a)', ' # Cells before Cgns_Mod_Merge_Nodes_New function (sample)'
    do c = 1, min(6, grid % n_cells)
      print '(a,8i8)', ' # ', &
        (grid % cells_n(i,c), i = 1, grid % cells_n_nodes(c))
    end do
  end if

  ! For each unique interface
  do int = 1, cnt_int

    ! Array new_seq is a map for nodes 'n -> new_seq(n)'
    allocate(new_seq (grid % n_nodes)); new_seq = 0
    allocate(nodes_to_remove(grid % n_nodes));  nodes_to_remove = .false.

    do n = 1, grid % n_nodes
      new_seq(n) = n
    end do ! n
    !--------------------------------------------------------------------------!
    !   new_seq map now is:                                                    !
    !   1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16                  !
    !--------------------------------------------------------------------------!

    if (verbose) then

      print '(a,i8)', ' # Interface to keep: ', int
      v = 0
      do c = 1, cnt_int_cells

        do k = 1, 4
          n1 = interface_cells(1, c, k, int)

          if (n1 > 0 .and. v < 7) then
            print '(a,4i8)', ' # ', interface_cells(1, c, 1:4, int)
            v = v + 1
          end if ! n1 > 0

        end do ! k
      end do ! c

      print '(a,i8)', ' # Interface to remove: ', int
      v = 0
      do c = 1, cnt_int_cells

        do k = 1, 4
          n2 = interface_cells(2, c, k, int)

          if (n2 > 0 .and. v < 7) then
            print '(a,4i8)', ' # ', interface_cells(2, c, 1:4, int)
            v = v + 1
          end if ! n1 > 0

        end do ! k
      end do ! c

    end if ! verbose

    ! Count nodes on interface
    cnt_nodes_on_int_to_keep = 0
    cnt_nodes_on_int_to_rem  = 0

    do c = 1, cnt_int_cells
      do k = 1, 4
        n1 = interface_cells(1, c, k, int)
        n2 = interface_cells(2, c, k, int)

        if (n1 > 0) cnt_nodes_on_int_to_keep = cnt_nodes_on_int_to_keep + 1
        if (n2 > 0) cnt_nodes_on_int_to_rem  = cnt_nodes_on_int_to_rem  + 1

      end do ! k
    end do ! c

    ! Allocate memory
    cnt_nodes_on_int_total = cnt_nodes_on_int_to_rem + cnt_nodes_on_int_to_keep
    allocate(criterion (cnt_nodes_on_int_total)); criterion = 0.
    allocate(old_seq   (cnt_nodes_on_int_total)); old_seq = 0

    !--------------------------------------!
    !   Prescribe some sorting criterion   !
    !--------------------------------------!
    cnt_node = 1 ! count grouped nodes to keep

    i = 1
    do c = 1, cnt_int_cells
      do k = 1, 4
        n = 0
        n1 = interface_cells(1, c, k, int)
        n2 = interface_cells(2, c, k, int)
        if (n1 > 0) n = n1
        if (n2 > 0) n = n2
        if (n > 0) then
          old_seq(i) = n
          criterion(i) = grid % xn(n) + grid % yn(n)*big + grid % zn(n)*big**2
          i = i + 1
        end if
      end do ! k
    end do ! n

    ! Sort nodes by this criterion
    call Sort_Mod_Real_Carry_Int(criterion, old_seq)

    ! Count grouped nodes after sorting
    v = 1 ! related to verbose output
    do i = 2, cnt_nodes_on_int_total
      ! If node is unique
      if( .not. Approx_Real(criterion(i-1), criterion(i), small) ) then
        cnt_node = cnt_node + 1
      end if
    end do ! i

    v = 1 ! related to verbose output
    do i = 2, cnt_nodes_on_int_total
      ! If node is duplicated
      if( Approx_Real(criterion(i-1), criterion(i), small) ) then

        if (old_seq(i-1) .ne. old_seq(i)) then
          new_seq(old_seq(i)  ) = minval(old_seq(i-1:i))
          new_seq(old_seq(i-1)) = minval(old_seq(i-1:i))

          nodes_to_remove(maxval(old_seq(i-1:i))) = .true.
        end if

        if (verbose .and. v <= min(6, grid % n_nodes)) then
          write (*, '(a)', advance='no')' # '
          print '(100a15)',('---------------', k = i-1, i)
          write (*, '(a)', advance='no')' # n: '
          print '(i13,a,i13)', minval(old_seq(i-1:i)), '<-', &
            maxval(old_seq(i-1:i))
          write (*, '(a)', advance='no')' # c: '
          print '(100es14.7)', (criterion(k), k = i-1, i)
          v = v + 1
        end if ! verbose
      end if
    end do ! i

    !--------------------------------------------------------------------------!
    !   new_seq map now is:                                                    !
    !   1  2  3  4  5  6  7  8  5  6  7  8  13  14  15  16                     !
    !--------------------------------------------------------------------------!

    !--------------------------------------------!
    !   Reconstruct new nodes and cells arrays   !
    !--------------------------------------------!
    cnt_node = 0
    do n = 1, grid % n_nodes
      if (nodes_to_remove(n)) cnt_node = cnt_node + 1
    end do
    cnt_node = grid % n_nodes - cnt_node

    allocate(x_new(cnt_node))
    allocate(y_new(cnt_node))
    allocate(z_new(cnt_node))

    cnt_node = 1
    do n = 1, grid % n_nodes
      if (.not. nodes_to_remove(n)) then ! if node is unique
        x_new(cnt_node) = grid % xn(n)
        y_new(cnt_node) = grid % yn(n)
        z_new(cnt_node) = grid % zn(n)

        cnt_node = cnt_node + 1
      else
        ! New sequence: shift all non-unique nodes in increasing order
        do c = cnt_node, grid % n_nodes
          if (new_seq(c) > cnt_node) then
            new_seq(c) = new_seq(c) -1
          end if
        end do ! c
      end if
    end do ! n
    !--------------------------------------------------------------------------!
    !   new_seq map now is:                                                    !
    !   1  2  3  4  5  6  7  8  5  6  7  8  9  10  11  12                      !
    !--------------------------------------------------------------------------!

    ! Remap nodes at cells according to new_seq
    do c = 1, grid % n_cells
      do i = 1, grid % cells_n_nodes(c)
        grid % cells_n(i,c) = new_seq(grid % cells_n(i,c))
      end do ! i
    end do ! c

    deallocate(criterion)
    deallocate(old_seq  )

    !------------------------!
    !   Reinitialize nodes   !
    !------------------------!
    print '(a,i7)', ' # New number of nodes: ', cnt_node - 1
    print '(a)', ' #-----------------------------------------------------------'

    deallocate(grid % xn)
    deallocate(grid % yn)
    deallocate(grid % zn)

    call Grid_Mod_Allocate_Nodes(grid, cnt_node - 1) ! -> grid % n_nodes

    grid % xn(1: grid % n_nodes) = x_new(1: grid % n_nodes)
    grid % yn(1: grid % n_nodes) = y_new(1: grid % n_nodes)
    grid % zn(1: grid % n_nodes) = z_new(1: grid % n_nodes)

    deallocate(x_new)
    deallocate(y_new)
    deallocate(z_new)

    !----------------------!
    !   Remap interfaces   !
    !----------------------!

    ! If not last interface, remap interfaces as well
    if (int < cnt_int) then

      do n = int + 1, cnt_int
        do c = 1, cnt_int_cells
          do k = 1, 4
            n1 = interface_cells(1, c, k, n)
            n2 = interface_cells(2, c, k, n)
            if (n1 > 0) interface_cells(1, c, k, n) = new_seq(n1)
            if (n2 > 0) interface_cells(2, c, k, n) = new_seq(n2)
          end do ! k
        end do ! c
      end do ! n

    end if

    deallocate(new_seq)
    deallocate(nodes_to_remove)

  end do ! interfaces

  if (verbose) then
    print *, '# Cells after Cgns_Mod_Merge_Nodes_New function (sample)'
    do c = 1, min(6, grid % n_cells)
      print '(a,8i8)', ' # ', &
        (grid % cells_n(i,c), i = 1, grid % cells_n_nodes(c))
    end do
  end if

  end subroutine
