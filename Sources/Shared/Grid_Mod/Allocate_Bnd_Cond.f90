!==============================================================================!
  subroutine Grid_Mod_Allocate_Bnd_Cond(grid)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!==============================================================================!

  ! Boundary conditions' data
  allocate(grid % bnd_cond % name(grid % n_bnd_cond))
  allocate(grid % bnd_cond % type(grid % n_bnd_cond))

  allocate(grid % bnd_cond % color_f(grid % n_bnd_cond))
  allocate(grid % bnd_cond % color_l(grid % n_bnd_cond))
  allocate(grid % bnd_cond % color(-grid % n_bnd_cells-1:-1))

  end subroutine
