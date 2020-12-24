!==============================================================================!
  subroutine Grid_Mod_Merge_Duplicate_Nodes(grid)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, e, i_nod, n, s, cnt
  integer, allocatable :: old_n(:)                ! node number
  integer, allocatable :: new_n(:)
  real,    allocatable :: xn(:), yn(:), zn(:)
!==============================================================================!

  print '(a)',    ' #================================================'
  print '(a)',    ' # Merging nodes if needed'
  print '(a)',    ' #------------------------------------------------'
  print '(a,i9)', ' # Original number of nodes:', grid % n_nodes

  !-----------------------------------------!
  !   Allocate memory for local variables   !
  !-----------------------------------------!
  allocate(xn   (grid % n_nodes))
  allocate(yn   (grid % n_nodes))
  allocate(zn   (grid % n_nodes))
  allocate(old_n(grid % n_nodes))
  allocate(new_n(grid % n_nodes))

  !-------------------------------!
  !   Sort nodes by coordinates   !
  !-------------------------------!

  ! Copy node coordinates
  do n = 1, grid % n_nodes
    old_n(n) = n
    xn(n) = grid % xn(n)
    yn(n) = grid % yn(n)
    zn(n) = grid % zn(n)
  end do

  ! Sort the nodes by their coordinates
  call Sort_Mod_3_Real_Carry_Int(xn, yn, zn, old_n)

  !-----------------------------------------------------!
  !   Look for duplicates and assign new node numbers   !
  !-----------------------------------------------------!
  cnt = 1
  new_n(old_n(1)) = cnt
  do n = 2, grid % n_nodes
    if( .not. (Math_Mod_Approx_Real(xn(n), xn(n-1)) .and.  &
               Math_Mod_Approx_Real(yn(n), yn(n-1)) .and.  &
               Math_Mod_Approx_Real(zn(n), zn(n-1))) ) then
      cnt = cnt + 1
    end if
    new_n(old_n(n)) = cnt
  end do
  print '(a,i9)', ' # Number of unique nodes:  ', cnt

  ! Decide what to do based on the compressed number of nodes
  if(cnt .eq. grid % n_nodes) then
    print '(a)', ' # No duplicate nodes found, nothing to merge!'
    return

  !---------------------------------!
  !   Do the actuall sorting work   !
  !---------------------------------!
  else
    print '(a,i0.0,a)', ' # ', grid % n_nodes - cnt, ' duplicate nodes' //  &
                        ' found; compressing them now!'

    call Sort_Mod_Real_By_Index(grid % n_nodes, grid % xn, new_n)
    call Sort_Mod_Real_By_Index(grid % n_nodes, grid % yn, new_n)
    call Sort_Mod_Real_By_Index(grid % n_nodes, grid % zn, new_n)

    ! Cells' nodes
    do c = -grid % n_bnd_cells, grid % n_cells
      do i_nod = 1, abs(grid % cells_n_nodes(c))
        n = grid % cells_n(i_nod, c)
        if(n > 0) grid % cells_n(i_nod, c) = new_n(n)
      end do
    end do

    ! Faces' nodes
    do s = 1, grid % n_faces
      do i_nod = 1, grid % faces_n_nodes(c)
        n = grid % faces_n(i_nod, s)
        if(n > 0) grid % faces_n(i_nod, s) = new_n(n)
      end do
    end do

    ! Edges' nodes
    do e = 1, grid % n_edges
      do i_nod = 1, 2
        n = grid % edges_n(i_nod, e)
        if(n > 0) grid % edges_n(i_nod, e) = new_n(n)
      end do
    end do

    print '(a)', ' # Done!'

  end if

  !-----------------------------------------!
  !   Deallocate memory used by variables   !
  !-----------------------------------------!
  deallocate(xn)
  deallocate(yn)
  deallocate(zn)
  deallocate(old_n)
  deallocate(new_n)

  end subroutine
