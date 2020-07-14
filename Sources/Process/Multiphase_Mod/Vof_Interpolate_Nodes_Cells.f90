!==============================================================================!
  subroutine Multiphase_Mod_Vof_Interpolate_Nodes_Cells(grid,                 &
                                                        var_node, var_cell)
!------------------------------------------------------------------------------!
!   Interpolates a variable from nodes to cell center                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  real            :: var_node(1 : grid % n_nodes)
  real            :: var_cell(-grid % n_bnd_cells:grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  integer :: n, s, c, c1, c2, i_nod, i_fac
  real    :: sum1, sum2
!==============================================================================!

  do c = 1, grid % n_cells
    sum1 = 0.0; sum2 = 0.0
    ! Loop on cells
    do i_nod = 1, grid % cells_n_nodes(c)
      n = grid % cells_n(i_nod, c)
      sum1 = sum1 + grid % cells_weight_n(i_nod, c) * var_node(n)
      sum2 = sum2 + grid % cells_weight_n(i_nod, c)
    end do
    var_cell(c) = sum1 / sum2
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, var_cell)

  end subroutine
