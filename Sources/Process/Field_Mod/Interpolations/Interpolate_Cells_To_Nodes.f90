!==============================================================================!
  subroutine Interpolate_Cells_To_Nodes(Flow, var_cell, var_node)
!------------------------------------------------------------------------------!
!   Interpolates a variable from cell centers to nodes                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type) :: Flow
  real, intent(in)  :: var_cell( -Flow % pnt_grid % n_bnd_cells  &
                                 :Flow % pnt_grid % n_cells)
  real, intent(out) :: var_node(1:Flow % pnt_grid % n_nodes)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: n, c, i_cel
!==============================================================================!

  ! Take alias
  Grid => Flow % pnt_grid

  ! It might need a call to: call Grid_Mod_Exchange_Cells_Real(Grid, var_cell)

  ! Interpolate from all cells to nodes
  do n = 1, Grid % n_nodes

    var_node(n) = 0.0

    ! Loop on cells
    do i_cel = 1, Grid % nodes_n_cells(n)
      c = Grid % nodes_c(i_cel, n)
      var_node(n) = var_node(n) + Grid % weight_c2n(i_cel, n) * var_cell(c)
    end do

  end do

  end subroutine
