!==============================================================================!
  integer function N_Nodes_In_Bnd_Color(Grid, bc, node_data)
!------------------------------------------------------------------------------!
!   Counts and marks (with node_data) nodes in the given boundary color        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: Grid
  integer         :: bc
  integer         :: node_data(Grid % n_nodes)
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, i_nod, n, cnt
!==============================================================================!

  ! Nullify on entry
  cnt = 0
  node_data(:) = 0

  ! Browse through all boundary cells and mark nodes
  ! of those cells in the given boundary condition
  do c = -Grid % n_bnd_cells, -1
    if( Grid % bnd_cond % color(c) .eq. bc ) then
      do i_nod = 1, Grid % cells_n_nodes(c)
        n = Grid % cells_n(i_nod, c)
        if(node_data(n) .eq. 0) then  ! hasn't been marked yet
          cnt = cnt + 1
          node_data(n) = cnt
        end if
      end do
    end if
  end do

  ! Return a sum of all marked nodes
  N_Nodes_In_Bnd_Color = cnt

  end function

