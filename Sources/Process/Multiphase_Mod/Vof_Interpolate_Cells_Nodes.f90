!==============================================================================!
  subroutine Multiphase_Mod_Vof_Interpolate_Cells_Nodes(grid,                 &
                                                        var_cell, var_node)
!------------------------------------------------------------------------------!
!   Interpolates a variable from cell centers to nodes                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  real            :: var_cell(-grid % n_bnd_cells:grid % n_cells)
  real            :: var_node(1:grid % n_nodes)
!-----------------------------------[Locals]-----------------------------------!
  integer :: n, s, c, c1, c2, i_cell, i_fac
  real    :: sum1, sum2
!==============================================================================!

  ! Interpolate from all cells to nodes
  do n = 1, grid % n_nodes
    sum1 = 0.0; sum2 = 0.0
    ! Loop on cells
    do i_cell = 1, grid % nodes_n_cells(n)
      c = grid % nodes_c(i_cell, n)
      sum1 = sum1 + grid % nodes_weight_c(i_cell, n) * var_cell(c)
      sum2 = sum2 + grid % nodes_weight_c(i_cell, n)
    end do
    var_node(n) = sum1 / sum2
  end do

  call Grid_Mod_Exchange_Nodes_Real(grid, var_node)

  end subroutine
