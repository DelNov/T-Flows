!==============================================================================!
  subroutine Interpolate_Nodes_To_Cells(Flow, var_node, var_cell)
!------------------------------------------------------------------------------!
!   Interpolates a variable from nodes to cells                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type)  :: Flow
  real, intent(in)  :: var_node(1:Flow % pnt_grid % n_nodes)
  real, intent(out) :: var_cell( -Flow % pnt_grid % n_bnd_cells  &
                                 :Flow % pnt_grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: n, c, i_nod
!==============================================================================!

  ! Take alias
  Grid => Flow % pnt_grid

  ! Interpolate from all nodes to cells
  do c = 1, Grid % n_cells

    var_cell(c) = 0.0

    ! Loop on nodes
    do i_nod = 1, abs(Grid % cells_n_nodes(c))

      n = Grid % cells_n(i_nod, c)
      var_cell(c) = var_cell(c)  &
                  + Grid % weight_n2c(i_nod, c) * var_node(n)

    end do

  end do

  call Grid_Mod_Exchange_Cells_Real(Grid, var_cell)

  end subroutine
