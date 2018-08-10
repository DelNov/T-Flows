!==============================================================================!
  subroutine Cgns_Mod_Merge_Nodes_Newest(grid)
!------------------------------------------------------------------------------!
!   Merges blocks by merging common surface without changing structure of      !
!   node and cells connection  tables                                          !
!                                                                              !
!   To do: put this procedure after reading block                              !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
  use Const_Mod, only: HUGE, PI
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!------------------------------------------------------------------------------!
  include "../Shared/Approx.int"
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, n, i, m, k, v, cnt_node, mn, mx
  
  real,    allocatable :: criterion(:) ! sorting criterion
  integer, allocatable :: old_seq(:), new_seq(:)
  integer, allocatable :: nodes_to_sift(:)
  real,    allocatable :: x_new(:), y_new(:), z_new(:)
  logical, allocatable :: nodes_to_remove(:) ! marked duplicated nodes to remove
  real                 :: big, small

  integer              :: int, cnt_nodes_to_remove
!==============================================================================!

  print *, '# Merging blocks since they have duplicating nodes '
  print *, '# Hint: Join blocks in mesh builder to avoid any problems'
  print '(a38,i9)', ' # Old number of nodes:               ', grid % n_nodes

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
  !   10:  0  1  0        cell 2                    cell 1                     !
  !   11:  0 -1  1    y                                                        !
  !   12:  0 -1  0    ^                                                        !
  !   13: -1  1  1    |   ^ z                                                  !
  !   14: -1  1  0    |  /                                                     !
  !   15: -1 -1  1    | /                                                      !
  !   16: -1 -1  0    --------> x                                              !
  !                                                                            !
  !   cell 1: 1   2   3   4   5   6   7   8                                    !
  !   cell 2: 9  10  11  12  13  14  15  16                                    !
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
  !   10: -1  1  0    cell 2                    cell 1                         !
  !   11: -1 -1  1                                                             !
  !   12: -1 -1  0                                                             !
  !                                                                            !
  !   cell 1: 1  2  3  4   5   6   7   8                                       !
  !   cell 2: 5  6  7  8   9  10  11  12                                       !
  !                                                                            !
  !----------------------------------------------------------------------------!

  print *, 'cells with dups:'
  do int = 1, cnt_int
    do c = 1, grid % n_cells
      if (.not. interface_cells(c,int)) then
        print *, 'int=', int, c
      end if
    end do
  end do

  print *, 'nodes to delete:'
  do int = 1, cnt_int
    do n = 1, grid % n_nodes
      if (.not. interface_nodes(n,int)) then
        print *, 'int=', int, n
      end if
    end do
  end do

  ! Count nodes to delete
  i = 1
  do int = 1, cnt_int ! for each unique interface

    do n = 1, grid % n_nodes
      if (interface_nodes(n,int)) then ! for each duplicated node
        i = i + 1
      end if
    end do

  end do

  cnt_nodes_to_remove = i - 1
  print *, 'nodes to delete:', cnt_nodes_to_remove

  !--------------------------------------------!
  !   Reconstruct new nodes and cells arrays   !
  !--------------------------------------------!
  cnt_node = grid % n_nodes - cnt_nodes_to_remove

  allocate(x_new(cnt_node))
  allocate(y_new(cnt_node))
  allocate(z_new(cnt_node))

  cnt_node = 1
  do n = 1, grid % n_nodes
    if (.not. interface_nodes(n,int)) then ! if node is unique

      x_new(cnt_node) = grid % xn(n)
      y_new(cnt_node) = grid % yn(n)
      z_new(cnt_node) = grid % zn(n)

      cnt_node = cnt_node + 1
    else
      ! New sequence: shift all non-unique nodes in decreasing order
      do c = cnt_node, grid % n_nodes
        if (new_seq(c) > cnt_node) then
          new_seq(c) = new_seq(c) -1
        end if
      end do ! c
    end if
  end do ! n

  print '(a38,i9)', ' # New number of nodes:               ', cnt_node - 1
  stop

  !----------------------------------------------------------------------------!
  !   new_seq now became:                                                      !
  !   1  2  3  4  5  6  7  8  5  6  7  8  9  10  11  12                        !
  !----------------------------------------------------------------------------!

  ! Remap nodes in cells according to new_seq
  do c = 1, grid % n_cells
    do i = 1, grid % cells_n_nodes(c)
      grid % cells_n(i,c) = new_seq(grid % cells_n(i,c))
    end do ! i
  end do ! c

  if (verbose) then
    print *, '# Cells after Cgns_Mod_Merge_Nodes_New function (sample)'
    do c = 1, 6
      print *, '#', (grid % cells_n(i,c), i = 1, grid % cells_n_nodes(c))
    end do
  end if

  cnt_node = cnt_node - 1

  !-----------------------!
  !   Reinitialize nodes  !
  !-----------------------!
  grid % n_nodes = cnt_node

  deallocate(grid % xn)
  deallocate(grid % yn)
  deallocate(grid % zn)

  call Grid_Mod_Allocate_Nodes(grid, grid % n_nodes)

  grid % xn(1: grid % n_nodes) = x_new(1: grid % n_nodes)
  grid % yn(1: grid % n_nodes) = y_new(1: grid % n_nodes)
  grid % zn(1: grid % n_nodes) = z_new(1: grid % n_nodes)

  deallocate(x_new)
  deallocate(y_new)
  deallocate(z_new)

  deallocate(criterion)
  deallocate(old_seq)
  deallocate(new_seq)
  deallocate(nodes_to_remove)

  end subroutine
