!==============================================================================!
  subroutine Cgns_Mod_Merge_Nodes_Old(grid)
!------------------------------------------------------------------------------!
!   Merges nodes from different blocks.                                        !
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
  integer              :: c, n, cnt_node
  real,    allocatable :: criter(:)          ! sorting criteria
  real,    allocatable :: x_old(:), y_old(:), z_old(:)
  integer, allocatable :: old_number(:)
  integer, allocatable :: new_number(:)
  integer, allocatable :: compressed(:)
  real                 :: big, small
!==============================================================================!

  print *, '# Old number of nodes: ', grid % n_nodes

  ! Allocate memory
  allocate(criter    (grid % n_nodes));  criter     = 0
  allocate(old_number(grid % n_nodes));  old_number = 0
  allocate(new_number(grid % n_nodes));  new_number = 0
  allocate(compressed(grid % n_nodes));  compressed = 0
  allocate(x_old(grid % n_nodes))
  allocate(y_old(grid % n_nodes))
  allocate(z_old(grid % n_nodes))
  
  ! Estimate big and small
  call Grid_Mod_Estimate_Big_And_Small(grid, big, small)

  !---------------------------!
  !   Store old coordinates   !
  !---------------------------!
  do n = 1, grid % n_nodes
    x_old(n) = grid % xn(n)
    y_old(n) = grid % yn(n)
    z_old(n) = grid % zn(n)
  end do

  !--------------------------------------!
  !   Prescribe some sorting criterion   !
  !--------------------------------------!
  do n = 1, grid % n_nodes
    criter(n) = grid % xn(n) + grid % yn(n) * big + grid % zn(n) * big * big
    old_number(n) = n
  end do

  !----------------------------------------------------------------------------!
  !   Original block structure with duplicate nodes:                           !
  !                                                                            !
  !   9--19--11--12 21--22--23--24                                             !
  !   |           |  |           |                                             !
  !   5   6   7   8 17  18  19  20                                             !
  !   |           |  |           |                                             !
  !   1---2---3---4 13--14--15--16                                             !
  !                                                                            !
  !   old_number:                                                              !
  !   1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24   !
  !----------------------------------------------------------------------------!

  !----------------------------------!
  !   Sort nodes by give criterion   !
  !----------------------------------!
  print *, '# Sorting the nodes ...'
  call Sort_Real_Carry_Int(criter(1), old_number(1), grid % n_nodes, 2)
  print *, '# ... done!'

  !----------------------------------------------------------------------------!
  !   After sorting by given criterion:                                        !
  !                                                                            ! 
  !   3---6---9--14 15--18--21--24                                             !
  !   |           |  |           |                                             !
  !   2   5   8  12 13  17  20  23                                             !
  !   |           |  |           |                                             !
  !   1---4---7--10 11--16--19--22                                             !
  !                                                                            !
  !   old_number:                                                              !
  !   1  5  9  2  6 19  3  7 11  4-13  8-17 12-21 14 18 22 15 19 23 16 20 24   !
  !----------------------------------------------------------------------------!

  !------------------------!
  !   Compress the nodes   !
  !------------------------!
  cnt_node = 1
  compressed(1) = cnt_node
  do n = 2, grid % n_nodes
    if( .not. Approx(criter(n), criter(n-1), SMALL) ) cnt_node = cnt_node + 1
    compressed(n) = cnt_node
  end do
  print *, '# New number of nodes: ', cnt_node

  !---------------------------------------------------------------------------!
  !   After compression:                                                      !
  !                                                                           !
  !   3---6---9--12--15--18--21                                               !
  !   |           |           |                                               !
  !   2   5   8  11  14  17  20                                               ! 
  !   |           |           |                                               !
  !   1---4---7--10--13--16--19                                               !
  !                                                                           !
  !   old_number:                                                             !
  !   1  5  9  2  6 19  3  7 11  4-13  8-17 12-21 14 18 22 15 19 23 16 20 24  !
  !                                                                           !
  !   compressed:                                                             ! 
  !   1  2  3  4  5  6  7  8  9 10-10 11-11 12-12 13 14 25 16 17 28 19 20 21  !
  !---------------------------------------------------------------------------!

  !--------------------------!
  !   Work out new numbers   !
  !--------------------------!
  do n = 1, grid % n_nodes
    new_number( old_number(n) ) = compressed(n)
  end do

  !-------------------------------------------------!
  !   Apply new node numbers to cells_n structure   !
  !-------------------------------------------------!
  do c = 1, grid % n_cells
    do n = 1, grid % cells_n_nodes(c)
      grid % cells_n(n,c) = new_number( grid % cells_n(n,c) )
    end do
  end do

  do n = 1, grid % n_nodes
    grid % xn( new_number(n) ) = x_old(n)
    grid % yn( new_number(n) ) = y_old(n)
    grid % zn( new_number(n) ) = z_old(n)
  end do

  !----------------------------------------------!
  !   Final touch: correct the number of nodes   !
  !----------------------------------------------!
  grid % n_nodes = cnt_node

  end subroutine
