!==============================================================================!
  subroutine Find_Edges(Grid)
!------------------------------------------------------------------------------!
!>  Finds all unique grid edges from face-node connectivity.  Each face
!>  contributes one edge between every pair of consecutive face nodes.  Edges are
!>  stored with increasing node numbers, sorted, and compressed to unique pairs.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: s, e, i_nod, j_nod, n1, n2
  integer              :: n_full_edges
  integer, allocatable :: full_edge_n(:,:)
!==============================================================================!

  if(allocated(Grid % edges_n))  deallocate(Grid % edges_n)
  if(allocated(Grid % edges_bc)) deallocate(Grid % edges_bc)
  if(allocated(Grid % edges_fb)) deallocate(Grid % edges_fb)

  !------------------------------------!
  !   Count all face-perimeter edges   !
  !------------------------------------!
  n_full_edges = 0
  do s = 1, Grid % n_faces
    n_full_edges = n_full_edges + Grid % faces_n_nodes(s)
  end do

  Assert(n_full_edges .gt. 0)

  allocate(full_edge_n(2, n_full_edges));  full_edge_n(:,:) = 0

  !-----------------------!
  !   Pile up all edges   !
  !-----------------------!
  e = 0
  do s = 1, Grid % n_faces
    do i_nod = 1, Grid % faces_n_nodes(s)

      j_nod = i_nod + 1
      if(j_nod > Grid % faces_n_nodes(s)) j_nod = 1

      n1 = Grid % faces_n(i_nod, s)
      n2 = Grid % faces_n(j_nod, s)

      Assert(n1 .gt. 0)
      Assert(n2 .gt. 0)
      Assert(n1 .ne. n2)

      e = e + 1
      full_edge_n(1, e) = min(n1, n2)
      full_edge_n(2, e) = max(n1, n2)
    end do
  end do

  Assert(e .eq. n_full_edges)

  !--------------------------------------!
  !   Sort by first and second node id   !
  !--------------------------------------!
  call Sort % Two_Int(full_edge_n(1, :), full_edge_n(2, :))

  !------------------------!
  !   Count unique edges   !
  !------------------------!
  Grid % n_edges = 1
  do e = 2, n_full_edges
    if(full_edge_n(1, e) .ne. full_edge_n(1, e-1) .or.  &
       full_edge_n(2, e) .ne. full_edge_n(2, e-1)) then
      Grid % n_edges = Grid % n_edges + 1
    end if
  end do

  !------------------------!
  !   Store unique edges   !
  !------------------------!
  allocate(Grid % edges_n(2, Grid % n_edges));  Grid % edges_n(:,:) = 0

  Grid % n_edges = 1
  Grid % edges_n(1, Grid % n_edges) = full_edge_n(1, 1)
  Grid % edges_n(2, Grid % n_edges) = full_edge_n(2, 1)

  do e = 2, n_full_edges
    if(full_edge_n(1, e) .ne. full_edge_n(1, e-1) .or.  &
       full_edge_n(2, e) .ne. full_edge_n(2, e-1)) then

      Grid % n_edges = Grid % n_edges + 1
      Grid % edges_n(1, Grid % n_edges) = full_edge_n(1, e)
      Grid % edges_n(2, Grid % n_edges) = full_edge_n(2, e)
    end if
  end do

  end subroutine
