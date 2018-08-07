!==============================================================================!
  subroutine Cgns_Mod_Merge_Nodes_New(grid)
!------------------------------------------------------------------------------!
!  Merges blocks by merging common surface without changing structure of node  !
!  and cells connection  tables                                                !
!  To do: put this procedure after reading block                               !
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!------------------------------------------------------------------------------!
  include "../Shared/Approx.int"
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, n, i, m, k, v, cnt_node, min, max
  real,    allocatable :: criterion(:) ! sorting criterion
  integer, allocatable :: old_seq(:), new_seq(:)
  real,    allocatable :: x_new(:), y_new(:), z_new(:)
  logical, allocatable :: nodes_to_remove(:) ! marked duplicated nodes to remove
  real,    parameter   :: BIG   = 104651.    ! prime number
  real,    parameter   :: SMALL = 1.e-6      ! precision for sorted criterion
!==============================================================================!

  print *, '# Merging blocks since they have duplicating nodes '
  print *, '# Hint: Join blocks in mesh builder to avoid any problems'
  print *, '# Old number of nodes: ', grid % n_nodes

  ! Allocate memory
  allocate(criterion      (grid % n_nodes));  criterion       = 0.
  allocate(old_seq        (grid % n_nodes));  old_seq         = 0
  allocate(new_seq        (grid % n_nodes));  new_seq         = 0
  allocate(nodes_to_remove(grid % n_nodes));  nodes_to_remove = .false.

  !--------------------------------------!
  !   Prescribe some sorting criterion   !
  !--------------------------------------!
  do n = 1, grid % n_nodes
    criterion(n) = grid % xn(n) + grid % yn(n) * BIG + grid % zn(n) * BIG ** 2
    old_seq(n) = n
  end do
  new_seq(:) = old_seq(:)

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
  !   old_seq:                                                                 !
  !    1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16                   !
  !----------------------------------------------------------------------------!

  !----------------------------------!
  !   Sort nodes by this criterion   !
  !----------------------------------!
  call Sort_Real_Carry_Int_Heapsort(criterion(1), old_seq(1), grid % n_nodes)

  if (verbose) then
    print *, '# Cells before Cgns_Mod_Merge_Nodes_New function (sample)'
    do c = 1, 6
      print *, '#', (grid % cells_n(i,c), i = 1, grid % cells_n_nodes(c))
    end do
  end if

  !----------------------------------------------------------------------------!
  !   old_seq is now:                                                          !
  !   16 12--8  4 14  6--10 2 15   7--11   3  13   5--9   1                    !
  !   "--" means that criterion has same value for these elements              !
  !----------------------------------------------------------------------------!
  cnt_node = 1
  if (verbose) print *,'# Uniting nodes (sample):'

  n = 1
  v = 1 ! related to verbose output

  do while ( n < grid % n_nodes )
    m = n + 1

    ! If node is unique
    if( .not. Approx(criterion(n), criterion(n+1), SMALL) ) then
      cnt_node = cnt_node + 1
    else ! if node is duplicated

      ! Check next nodes on the list by criterion
      do while (Approx(criterion(m), criterion(m+1), SMALL))
        m = m + 1
      end do
      ! [n : m] are duplicated

      min = minval(old_seq(n:m))
      max = maxval(old_seq(n:m))

      ! Mark node to remove
      do k = n, m
        if ( old_seq(k) .ne. min ) then
          nodes_to_remove(old_seq(k)) = .true.
          ! New sequence (stage 1) : substitute duplicated modes unique
          new_seq(old_seq(k)) = min
        end if
      end do

      if (verbose .and. v < 7) then
        write (*, '(a)', advance='no')' # '
        print '(100a15)',('---------------', k = n, m)
        write (*, '(a)', advance='no')' # n: '
        print '(100i14)', (old_seq(k), k = n, m)
        write (*, '(a)', advance='no')' # c: '
        print '(100es14.7)', (criterion(k), k = n, m)
        v = v + 1
      end if

    end if
    n = m
  end do

  !----------------------------------------------------------------------------!
  !   new_seq now became:                                                      !
  !   1  2  3  4  5  6  7  8  5  6  7  8  13  14  15  16                       !
  !----------------------------------------------------------------------------!

  !--------------------------------------------!
  !   Reconstruct new nodes and cells arrays   !
  !--------------------------------------------!

  allocate(x_new(1:cnt_node))
  allocate(y_new(1:cnt_node))
  allocate(z_new(1:cnt_node))

  cnt_node = 1
  do n = 1, grid % n_nodes
    if (.not. nodes_to_remove(n)) then ! if node is unique

      x_new(cnt_node) = grid % xn(n)
      y_new(cnt_node) = grid % yn(n)
      z_new(cnt_node) = grid % zn(n)

      cnt_node = cnt_node + 1
    else
      ! New sequence (stage 2) : shift all non-unique nodes in increasing order
      do c = cnt_node, grid % n_nodes
        if (new_seq(c) > cnt_node) then
          new_seq(c) = new_seq(c) -1
        end if
      end do ! c
    end if
  end do ! n

  print *, '# New number of nodes: ', cnt_node - 1

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