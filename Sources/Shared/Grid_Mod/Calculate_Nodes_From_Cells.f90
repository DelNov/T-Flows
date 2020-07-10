!==============================================================================!
  subroutine Grid_Mod_Calculate_Nodes_From_Cells(grid, var_cell, var_node)
!------------------------------------------------------------------------------!
!   Calculates nodal from cell values.                                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  real            :: var_cell(-grid % n_bnd_cells : grid % n_cells)
  real            :: var_node( 1 : grid % n_nodes)
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, n, i_cell
!==============================================================================!

  !--------------------------------------------------!
  !   Browse through all nodes to calculate values   !
  !   in the nodes from cell values surrounding it   !
  !--------------------------------------------------!
  do n = 1, grid % n_nodes
    var_node(n) = 0.0
    do i_cell = 1, grid % nodes_n_cells(n)
      c = grid % nodes_c(i_cell, n)
      var_node(n) = var_node(n)  &
                  + var_cell(c) * grid % nodes_weight_c(i_cell, n)
    end do
    var_node(n) = var_node(n)  &
                / sum(grid % nodes_weight_c(1:grid % nodes_n_cells(n), n))
  end do

  end subroutine
