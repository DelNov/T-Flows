!==============================================================================!
  integer function N_Nodes_At_Boundary(Convert, Grid, node_flag_at_boundary)
!------------------------------------------------------------------------------!
!>  Counts and marks nodes associated with all boundarues.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert  !! parent class
  type(Grid_Type)     :: Grid     !! grid being converted
  integer             :: node_flag_at_boundary(Grid % n_nodes)
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, i_nod, n, cnt
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Convert)
!==============================================================================!

  ! Nullify on entry
  cnt = 0
  node_flag_at_boundary(:) = 0

  ! Browse through all boundary cells and mark nodes
  ! of those cells in the given boundary condition
  do c = -Grid % n_bnd_cells, -1
    do i_nod = 1, Grid % cells_n_nodes(c)
      n = Grid % cells_n(i_nod, c)
      if(node_flag_at_boundary(n) .eq. 0) then  ! hasn't been marked yet
        cnt = cnt + 1
        node_flag_at_boundary(n) = 1
      end if
    end do
  end do

  ! Return a sum of all marked nodes
  N_Nodes_At_Boundary = cnt

  end function

