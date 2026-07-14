!==============================================================================!
  integer function N_Sharp_Corners(Convert, Grid,  &
                                   sharp_edge_flag, sharp_node_rank)
!------------------------------------------------------------------------------!
!>  This function is designed to identify and count sharp corners in a 3D
!>  grid structure. It takes the Grid as an input and returns an array
!>  sharp_node_rank marking the nodes at sharp corners.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type)  :: Convert
  type(Grid_Type)      :: Grid
  integer, intent(in)  :: sharp_edge_flag(Grid % n_edges)
  integer, intent(out) :: sharp_node_rank(Grid % n_nodes)
!-----------------------------------[Locals]-----------------------------------!
  integer :: e, cnt, n1, n2, n
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Convert)
!==============================================================================!

  ! Nullify on entry
  sharp_node_rank(:) = 0

  !------------------------------------------!
  !   Assign marks to nodes at sharp edges   !
  !------------------------------------------!
  do e = 1, Grid % n_edges

    ! If it is sharp, both convex (+1) and concave (-1)
    if(sharp_edge_flag(e) .ne. 0) then

      ! Take edges' nodes and mark them (increase
      ! the number of times they have been visited)
      n1 = Grid % edges_n(1, e)
      n2 = Grid % edges_n(2, e)

      sharp_node_rank(n1) = sharp_node_rank(n1) + 1
      sharp_node_rank(n2) = sharp_node_rank(n2) + 1
    end if

  end do

  !----------------------------------------------------------------------!
  !   Assign ranks to nodes which have been marked more than two times   !
  !----------------------------------------------------------------------!
  cnt = 0
  do n = 1, Grid % n_nodes
    if(sharp_node_rank(n) .gt. 2) then
      cnt = cnt + 1
      sharp_node_rank(n) = cnt
    else
      sharp_node_rank(n) = 0
    end if
  end do

  ! Return a sum of all marked nodes
  N_Sharp_Corners = cnt

  end function

