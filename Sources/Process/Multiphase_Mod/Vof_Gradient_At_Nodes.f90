!==============================================================================!
  subroutine Multiphase_Mod_Vof_Gradient_At_Nodes(grid, var_node, var_cell,   &
                                                  grad_x, grad_y, grad_z)
!------------------------------------------------------------------------------!
!   Computes gradient at nodes, based on cell values                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)       ,target :: grid
  real                          :: var_node   (1 : grid % n_nodes),   &
                                   var_cell   (-grid % n_bnd_cells    &
                                              : grid % n_cells),      &
                                   grad_x     (1 : grid % n_nodes),   &
                                   grad_y     (1 : grid % n_nodes),   &
                                   grad_z     (1 : grid % n_nodes)
!-----------------------------------[Locals]-----------------------------------!
  real,                 pointer :: wx(:,:), wy(:,:), wz(:,:)
  integer                       :: c, n, i_cell
!==============================================================================!

  wx => grid % weight_gradx_cells
  wy => grid % weight_grady_cells
  wz => grid % weight_gradz_cells

  grad_x = 0.0
  grad_y = 0.0
  grad_z = 0.0

  do n = 1, grid % n_nodes
    ! Loop on cells
    do i_cell = 1, grid % nodes_n_cells(n)
      c = grid % nodes_c(i_cell, n)
      grad_x(n) = grad_x(n) + wx(i_cell, n) * (var_cell(c) - var_node(n))
      grad_y(n) = grad_y(n) + wy(i_cell, n) * (var_cell(c) - var_node(n))
      grad_z(n) = grad_z(n) + wz(i_cell, n) * (var_cell(c) - var_node(n))
    end do
  end do

  call Grid_Mod_Exchange_Nodes_Real(grid, grad_x)
  call Grid_Mod_Exchange_Nodes_Real(grid, grad_y)
  call Grid_Mod_Exchange_Nodes_Real(grid, grad_z)

  end subroutine
