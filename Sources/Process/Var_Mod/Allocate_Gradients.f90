!==============================================================================!
  subroutine Var_Mod_Allocate_Gradients(phi)
!------------------------------------------------------------------------------!
!   This is to allocate additional values for statistics.                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type) :: phi
!-----------------------------------[Locals]-----------------------------------!
  integer :: n_cells, n_bnd_cells
!==============================================================================!

  ! Fetch resolutions
  n_bnd_cells = phi % pnt_grid % n_bnd_cells
  n_cells     = phi % pnt_grid % n_cells

  ! Gradients
  allocate (phi % x(-n_bnd_cells: n_cells));  phi % x = 0.
  allocate (phi % y(-n_bnd_cells: n_cells));  phi % y = 0.
  allocate (phi % z(-n_bnd_cells: n_cells));  phi % z = 0.

  end subroutine
