!==============================================================================!
  subroutine Var_Mod_Allocate_Statistics(phi)
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

  ! Terms for statistics
  allocate (phi % mean(-n_bnd_cells: n_cells));  phi % mean = 0.
  allocate (phi % fluc(-n_bnd_cells: n_cells));  phi % fluc = 0.
  allocate (phi % filt(-n_bnd_cells: n_cells));  phi % filt = 0.

  end subroutine
