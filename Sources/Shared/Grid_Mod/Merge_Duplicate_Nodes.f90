!==============================================================================!
  subroutine Merge_Duplicate_Nodes(Grid)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, e, i_nod, n, s, cnt
  real, allocatable :: xn(:), yn(:), zn(:)
!==============================================================================!

  if(First_Proc()) then
    print '(a)',     ' #======================================================='
    print '(a,a,a)', ' # Merging nodes in "', trim(Grid % name), '" if needed'
    print '(a)',     ' #-------------------------------------------------------'
    print '(a,i9)',  ' # Original number of nodes:', Grid % n_nodes
  end if

  !-----------------------------------------!
  !   Allocate memory for local variables   !
  !-----------------------------------------!
  allocate(xn(Grid % n_nodes))
  allocate(yn(Grid % n_nodes))
  allocate(zn(Grid % n_nodes))

  !-------------------------------!
  !   Sort nodes by coordinates   !
  !-------------------------------!

  ! Copy node coordinates
  do n = 1, Grid % n_nodes
    Grid % old_n(n) = n
    xn(n) = Grid % xn(n)
    yn(n) = Grid % yn(n)
    zn(n) = Grid % zn(n)
  end do

  ! Sort the nodes by their coordinates
  call Sort % Three_Real_Carry_Int(xn, yn, zn, Grid % old_n)

  !-----------------------------------------------------!
  !   Look for duplicates and assign new node numbers   !
  !-----------------------------------------------------!
  cnt = 1
  Grid % new_n(Grid % old_n(1)) = cnt
  do n = 2, Grid % n_nodes
    if( .not. (Math % Approx_Real(xn(n), xn(n-1)) .and.  &
               Math % Approx_Real(yn(n), yn(n-1)) .and.  &
               Math % Approx_Real(zn(n), zn(n-1))) ) then
      cnt = cnt + 1
    end if
    Grid % new_n(Grid % old_n(n)) = cnt
  end do
  if(First_Proc()) then
    print '(a,i9)', ' # Number of unique nodes:  ', cnt
  end if

  ! Decide what to do based on the compressed number of nodes
  if(cnt .eq. Grid % n_nodes) then
    if(First_Proc()) then
      print '(a)', ' # No duplicate nodes found, nothing to merge!'
    end if
    return

  !---------------------------------!
  !   Do the actuall sorting work   !
  !---------------------------------!
  else
    if(First_Proc()) then
      print '(a,i0.0,a)', ' # ', Grid % n_nodes - cnt, ' duplicate nodes' //  &
                          ' found; compressing them now!'
    end if

    call Sort % Real_By_Index(Grid % n_nodes, Grid % xn, Grid % new_n)
    call Sort % Real_By_Index(Grid % n_nodes, Grid % yn, Grid % new_n)
    call Sort % Real_By_Index(Grid % n_nodes, Grid % zn, Grid % new_n)

    ! Cells' nodes
    do c = -Grid % n_bnd_cells, Grid % n_cells
      do i_nod = 1, abs(Grid % cells_n_nodes(c))
        n = Grid % cells_n(i_nod, c)
        if(n > 0) Grid % cells_n(i_nod, c) = Grid % new_n(n)
      end do
    end do

    ! Faces' nodes
    do s = 1, Grid % n_faces
      do i_nod = 1, Grid % faces_n_nodes(c)
        n = Grid % faces_n(i_nod, s)
        if(n > 0) Grid % faces_n(i_nod, s) = Grid % new_n(n)
      end do
    end do

    ! Edges' nodes
    do e = 1, Grid % n_edges
      do i_nod = 1, 2
        n = Grid % edges_n(i_nod, e)
        if(n > 0) Grid % edges_n(i_nod, e) = Grid % new_n(n)
      end do
    end do

    ! Change the field for number of nodes
    ! (How come I didn't do it before?)
    Grid % n_nodes = 0
    do c = -Grid % n_bnd_cells, Grid % n_cells     ! cells' nodes
      do i_nod = 1, abs(Grid % cells_n_nodes(c))
        n = Grid % cells_n(i_nod, c)
        Grid % n_nodes = max(Grid % n_nodes, n)
      end do
    end do
    do s = 1, Grid % n_faces                       ! faces' nodes
      do i_nod = 1, Grid % faces_n_nodes(c)
        n = Grid % faces_n(i_nod, s)
        Grid % n_nodes = max(Grid % n_nodes, n)
      end do
    end do
    do e = 1, Grid % n_edges                       ! edges' nodes
      do i_nod = 1, 2
        n = Grid % edges_n(i_nod, e)
        Grid % n_nodes = max(Grid % n_nodes, n)
      end do
    end do

    if(First_Proc()) then
      print '(a)', ' # Done!'
    end if

  end if

  end subroutine
