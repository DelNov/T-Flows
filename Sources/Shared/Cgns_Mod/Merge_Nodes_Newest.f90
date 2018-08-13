!==============================================================================!
  subroutine Cgns_Mod_Merge_Nodes_Newest(grid)
!------------------------------------------------------------------------------!
!   For each interface in geometry merges nodes on interfaces and remaps       !
!   cell_connections.                                                          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!------------------------------------------------------------------------------!
  include "../Shared/Approx.int"
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, cnt_node, n, i, k, v
  
  real,    allocatable :: criterion(:) ! sorting criterion
  integer, allocatable :: old_seq(:), new_seq(:), interface_tmp(:,:)
  real,    allocatable :: x_new(:), y_new(:), z_new(:)
  real                 :: big, small
  integer              :: int, cnt_nodes_on_interface
!==============================================================================!

  print *, '# Merging blocks since they have common interfaces'
  print *, '# Hint: Join blocks in mesh builder to avoid any problems'
  print '(a38,i9)', ' # Old number of nodes:               ', grid % n_nodes
  
  !----------------------------------------------------------------------------!
  !   At this point number of interfaces cnt_int is known.                     !
  !   Nodes on interfaces are marked by interface_nodes .ne. 0                 !
  !   Value -1 of interface_nodes means this node has to be removed as         !
  !   duplicated, cell connection - renumbered.                                !
  !   Unfortunately coordinates of two nodes on same interface                 !
  !   interface_nodes(n, int) = 1 and interface_nodes(n, int) = -1             !
  !   do not match, because they originated from different blocks.             !
  !   That is why sorting is required by some criterion to group them.         !
  !   Fortunately value of unique nodes on each interface is known and         !
  !   determined by interface_nodes .eq. 1 or .eq. -1 .                        !
  !----------------------------------------------------------------------------!
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

  ! Estimate big and small
  call Grid_Mod_Estimate_Big_And_Small(grid, big, small)

  if (verbose) then
    print *, '# Cells before Cgns_Mod_Merge_Nodes_New function (sample)'
    do c = 1, min(6, grid % n_cells)
      print *, '#', (grid % cells_n(i,c), i = 1, grid % cells_n_nodes(c))
    end do
  end if

  ! For each unique interface
  do int = 1, cnt_int
  
    ! Array new_seq is a map for nodes "n -> new_seq(n)"
    allocate(new_seq (grid % n_nodes)); new_seq = 0
    do n = 1, grid % n_nodes
      new_seq(n) = n
    end do
    !--------------------------------------------------------------------------!
    !   new_seq map now is:                                                    !
    !   1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16                  !
    !--------------------------------------------------------------------------!

    if (verbose) then
      i = 1
      k = 1
      do n = 1, grid % n_nodes
        if ( interface_nodes(n, int) .eq. -1 ) then
          i = i + 1
          !print *, 'int=', int, ' n = ', n, 'del=true'
        else if ( interface_nodes(n, int) .eq. 1 ) then
          k = k + 1
          !print *, 'int=', int, ' n = ', n, 'del=false'
        end if
      end do
      print *, '# Interface', int, ', Nodes to remove:', i - 1, &
        'Nodes to keep:', k - 1
    end if

    ! Count nodes on interface
    cnt_node = 1
    do n = 1, grid % n_nodes
      if (interface_nodes(n, int) .ne. 0) then
        cnt_node = cnt_node + 1
      end if
    end do
    cnt_nodes_on_interface = cnt_node - 1
    
    !----------------------------------------------!
    !   Sind nodes until target value is reached   !
    !----------------------------------------------!
    if (verbose) print *, "Sind until n:", cnt_nodes_on_interface/2
    
    do while ( cnt_node .ne. cnt_nodes_on_interface/2 )
      cnt_node = 1

      ! Allocate memory
      allocate(criterion (cnt_nodes_on_interface)); criterion = 0.
      allocate(old_seq   (cnt_nodes_on_interface)); old_seq = 0
      
      !--------------------------------------!
      !   Prescribe some sorting criterion   !
      !--------------------------------------!

      i = 1
      do n = 1, grid % n_nodes
        if (interface_nodes(n, int) .ne. 0) then
          old_seq(i) = n
          criterion(i) = grid % xn(n) + grid % yn(n)*big + grid % zn(n)*big**2
          i = i + 1
        end if
      end do
      
      ! Sort nodes by this criterion
      call Sort_Real_Carry_Int_Heapsort(criterion(1), old_seq(1), &
      cnt_nodes_on_interface)
      
      ! Count grouped nodes after sorting
      v = 1 ! related to verbose output
      do i = 2, cnt_nodes_on_interface
        ! If node is unique
        if( .not. Approx(criterion(i-1), criterion(i), small) ) then
          cnt_node = cnt_node + 1
        end if
      end do ! i

      if (verbose) print *, '# Interface', int, ', Reached n:', cnt_node
      if (cnt_node < cnt_nodes_on_interface/2) then
        small = small / 2
        deallocate(criterion)
        deallocate(old_seq  )
      elseif (cnt_node > cnt_nodes_on_interface/2) then
        small = small * 2
        deallocate(criterion)
        deallocate(old_seq  )
      end if

    end do ! sind
    ! Now target value of nodes on interface is reached

    v = 1 ! related to verbose output
    do i = 2, cnt_nodes_on_interface
      ! If node is duplicated
      if( Approx(criterion(i-1), criterion(i), small) ) then

        ! Duplicated nodes must be assinged their lowest node id
        !if (interface_nodes(old_seq(i-1), int) .eq. -1) then
        !  new_seq(old_seq(i-1)) = new_seq(old_seq(i))
        !else
        !  new_seq(old_seq(i)) = new_seq(old_seq(i-1))
        !end if

        new_seq(old_seq(i)  ) = minval(old_seq(i-1:i))
        new_seq(old_seq(i-1)) = minval(old_seq(i-1:i))

        if (verbose .and. v <= min(6, grid % n_nodes)) then
          write (*, '(a)', advance='no')' # '
          print '(100a15)',('---------------', k = i-1, i)
          write (*, '(a)', advance='no')' # n: '
          print '(i13,a,i13)', minval(old_seq(i-1:i)), '<-', &
            maxval(old_seq(i-1:i))
          write (*, '(a)', advance='no')' # c: '
          print '(100es14.7)', (criterion(k), k = i-1, i)
          v = v + 1
        end if
      end if
    end do ! i

    !--------------------------------------------------------------------------!
    !   new_seq map now is:                                                    !
    !   1  2  3  4  5  6  7  8  5  6  7  8  13  14  15  16                     !
    !--------------------------------------------------------------------------!

    !--------------------------------------------!
    !   Reconstruct new nodes and cells arrays   !
    !--------------------------------------------!
    cnt_node = grid % n_nodes - cnt_nodes_on_interface/2
  
    allocate(x_new(cnt_node))
    allocate(y_new(cnt_node))
    allocate(z_new(cnt_node))

    cnt_node = 1
    do n = 1, grid % n_nodes
      if (interface_nodes(n, int) > -1) then

        x_new(cnt_node) = grid % xn(n)
        y_new(cnt_node) = grid % yn(n)
        z_new(cnt_node) = grid % zn(n)

        cnt_node = cnt_node + 1
      else
        do c = cnt_node, grid % n_nodes
          if (new_seq(c) > cnt_node) then
            new_seq(c) = new_seq(c) - 1
          end if
        end do ! c
      end if
    end do ! n

    ! If not last interface, remap interfaces as well
    if (int < cnt_int) then

      allocate(interface_tmp(cnt_node - 1, cnt_node)); interface_tmp = 0

      do k = int + 1, cnt_int

        cnt_node = 1
        do n = 1, grid % n_nodes
            if (interface_nodes(n, int) > -1) then
              interface_tmp(cnt_node, k) = interface_nodes(n, k)
              cnt_node = cnt_node + 1
            end if
        end do ! k
      end do ! n
    end if ! if not last interface

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
    print '(a38,i9)', ' # New number of nodes:               ', cnt_node - 1

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

    !-----------------------------!
    !   Reinitialize interfaces   !
    !-----------------------------!

    ! If not last interface, remap interfaces as well
    if (int < cnt_int) then
      deallocate(interface_nodes)
      allocate(interface_nodes(1:grid % n_nodes, cnt_int)); interface_nodes = 0

      do k = int + 1, cnt_int
        interface_nodes(1: grid % n_nodes, k) =  &
          interface_tmp(1: grid % n_nodes, k)
      end do ! k

      deallocate(interface_tmp)
    end if

    deallocate(new_seq)

  end do ! interfaces

  if (verbose) then
    print *, '# Cells after Cgns_Mod_Merge_Nodes_New function (sample)'
    do c = 1, min(6, grid % n_cells)
      print *, '#', (grid % cells_n(i,c), i = 1, grid % cells_n_nodes(c))
    end do
  end if

  end subroutine