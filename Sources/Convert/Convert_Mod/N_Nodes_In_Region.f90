!==============================================================================!
  integer function N_Nodes_In_Region(Convert, Grid, bc, node_rank_in_region)
!------------------------------------------------------------------------------!
!>  Counts and marks nodes associated with a specific boundary condition
!>  within a grid.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert  !! parent class
  type(Grid_Type)     :: Grid     !! grid being converted
  integer             :: bc       !! boundary condition rank (number)
  integer             :: node_rank_in_region(Grid % n_nodes)
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, i_nod, n, cnt
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Convert)
!==============================================================================!

  ! Nullify on entry
  cnt = 0
  node_rank_in_region(:) = 0

  ! Browse through all boundary cells and mark nodes
  ! of those cells in the given boundary condition
  do c = -Grid % n_bnd_cells, -1
    if( Grid % region % at_cell(c) .eq. bc ) then
      do i_nod = 1, Grid % cells_n_nodes(c)
        n = Grid % cells_n(i_nod, c)
        if(node_rank_in_region(n) .eq. 0) then  ! hasn't been marked yet
          cnt = cnt + 1
          node_rank_in_region(n) = cnt
        end if
      end do
    end if
  end do

  ! Return a sum of all marked nodes
  N_Nodes_In_Region = cnt

  end function

