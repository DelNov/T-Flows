!==============================================================================!
  subroutine Field_Mod_Interpolate_Nodes_To_Cells(flow, var_node, var_cell)
!------------------------------------------------------------------------------!
!   Interpolates a variable from nodes to cells                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type)  :: flow
  real, intent(in)  :: var_node(1:flow % pnt_grid % n_nodes)
  real, intent(out) :: var_cell( -flow % pnt_grid % n_bnd_cells  &
                                 :flow % pnt_grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: n, s, c, c1, c2, i_nod
!==============================================================================!

  ! Take alias
  grid => flow % pnt_grid

  do c = 1, grid % n_cells

    var_cell(c) = 0.0

    ! Loop on nodes
    do i_nod = 1, abs(grid % cells_n_nodes(c))

      n = grid % cells_n(i_nod, c)
      var_cell(c) = var_cell(c)  &
                  + grid % weight_n2c(i_nod, c) * var_node(n)

    end do

  end do

  call Grid_Mod_Exchange_Cells_Real(grid, var_cell)

  end subroutine
