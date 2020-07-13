!==============================================================================!
  subroutine Multiphase_Mod_Vof_Nodal_Gradient(grid, var_cell, var_node,   &
                                               grad_x, grad_y, grad_z)
!------------------------------------------------------------------------------!
!   Computes gradient at cell, based on nodal values                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)       ,target :: grid
  real                          :: var_cell   (-grid % n_bnd_cells    &
                                              : grid % n_cells),      &
                                   var_node   (1 : grid % n_nodes),   &
                                   grad_x     (-grid % n_bnd_cells    &
                                              : grid % n_cells),      &
                                   grad_y     (-grid % n_bnd_cells    &
                                              : grid % n_cells),      &
                                   grad_z     (-grid % n_bnd_cells    &
                                              : grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  real,                 pointer :: wx(:,:), wy(:,:), wz(:,:)
  integer                       :: c, n, i_nod
!==============================================================================!

  wx => grid % weight_gradx_nodes
  wy => grid % weight_grady_nodes
  wz => grid % weight_gradz_nodes

  grad_x = 0.0
  grad_y = 0.0
  grad_z = 0.0

  do c = 1, grid % n_cells
    ! Loop on nodes
    do i_nod = 1, grid % cells_n_nodes(c)
      n = grid % cells_n(i_nod, c)
      grad_x(c) = grad_x(c) + wx(i_nod, c) * (var_node(n) - var_cell(c))
      grad_y(c) = grad_y(c) + wy(i_nod, c) * (var_node(n) - var_cell(c))
      grad_z(c) = grad_z(c) + wz(i_nod, c) * (var_node(n) - var_cell(c))
    end do
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, grad_x)
  call Grid_Mod_Exchange_Cells_Real(grid, grad_y)
  call Grid_Mod_Exchange_Cells_Real(grid, grad_z)

  end subroutine
