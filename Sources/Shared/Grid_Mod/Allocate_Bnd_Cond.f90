!==============================================================================!
  subroutine Allocate_Bnd_Cond(Grid)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
!==============================================================================!

  ! Boundary conditions' data
  allocate(Grid % bnd_cond % name(Grid % n_bnd_cond))
  allocate(Grid % bnd_cond % type(Grid % n_bnd_cond))

  allocate(Grid % bnd_cond % color_f(Grid % n_bnd_cond))
  allocate(Grid % bnd_cond % color_l(Grid % n_bnd_cond))
  allocate(Grid % bnd_cond % color(-Grid % n_bnd_cells-1:-1))

  end subroutine
